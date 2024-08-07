function MODEL = fanchart(MODEL, varargin)

%{
    Genera las gráficas de las variables que son parte del Modelo. Pueden
    compararse, o no, con otro escenario o con otro corrimiento.
{
## Syntax ##

    MODEL = simPlots(MODEL, varargin)

## Input Arguments ##

__`MODEL`__ [ struct ] -
Debe contener al menos la estructura con los resultados del proceso de
simulación MODEL.F_pred.

* 'StartDate' = {} [ `Cell` ] - fechas de inicio del plot (pueden ser una o mas).

* 'EndDatePlot' = {} [ `Cell` ] - fechas de fin del plot (pueden ser una o mas).

* 'PlotList' = {'D4L_CPI'} [ `Cell` ] - Lista de variables a graficar. Si
no se le da la lista genera por default la de inflación total interanual


* 'apertura' = [] [ `Cell` ] - array con Vectores COLUMNA con valores de ajuste para el error
estándar (re-escala el error estándar por trimestre), por cada variable de la lista PlotList. 
Abre/cierra el abanico de forma Simétrica. Debe contener el mismo número 
de valores que trimestres de abanico.


* 'sesgo' = [] [ `Cell` ] - array con Vectores FILA con valores de ajuste para los valores de
los percentiles (re-escala la distribución), por cada variable de la lista PlotList.
Aumento de sesgo: valores mayores a 1.
Reducción del sesgo: valores menores a 1.

El vector debe ser del mismo tamaño que la cantidad de percentiles en la distribución.

En el caso del sesgo al alza/baja se ajustan los valores del vector que están
por ENCIMA/DEBAJO del percentil 50.
Por ejemplo, para perc del 10 al 90 en pasos de 0.10:

[1 1 1 1 1 1.25 1.25 1.25 1.25] para un fanchart con sesgo al alza (Asimetría hacia arriba).
[1.25 1.25 1.25 1.25 1 1 1 1 1] para un fanchart con sesgo a la baja (Asimetría hacia abajo).

En ambos casos pueden utilizarse factores distintos en cada elemento del
vector para ajustar la forma del abanico.


- DIE
- Agosto 2024
- MJGM
%}
p = inputParser;
    addParameter(p, 'StartDatePlot', MODEL.DATES.hist_end-20);
    addParameter(p, 'EndDatePlot', MODEL.DATES.pred_start+7);
    addParameter(p, 'SavePath', {});
    addParameter(p, 'Esc', {});
    addParameter(p, 'PlotList', {'D4L_CPI'});
    addParameter(p, 'Grilla', 0.10:0.10:0.90);
    addParameter(p, 'apertura', {});
    addParameter(p, 'sesgo', {});
    addParameter(p, 'CloseAll', true);    
    addParameter(p, 'AutoSave', true);
parse(p, varargin{:});
params = p.Results;

%SavePath
if  isempty(params.Esc)
    params.SavePath = fullfile('plots', MODEL.CORR_DATE, 'v0', 'fanchart');
else
    params.SavePath = fullfile('plots', MODEL.CORR_DATE, params.Esc, 'fanchart');
end

% Verificación y creación del directorio para las gráficas
if ~isfolder(params.SavePath)
    mkdir(params.SavePath)
else
    rmdir(params.SavePath, 's')
    mkdir(params.SavePath)
end  

%% Cálculo de Insumos

if isempty(params.apertura)
    for j = 1:length(params.PlotList)
        params.apertura{j} = {};
    end
end
if isempty(params.sesgo)
    for j = 1:length(params.PlotList)
        params.sesgo{j} = {};
    end
end


% Cálculo del Forecast Mean Squared Error
[~,~,FMSE] = fmse(MODEL.MF,...
                        length(MODEL.DATES.pred_start:params.EndDatePlot),...
                        'Select',params.PlotList);
FMSE = FMSE*params.PlotList;

% Fechas para cálculos del abanico (8 períodos)
MODEL.DATES.Fanchart   = MODEL.DATES.pred_start:params.EndDatePlot;
% Fechas para gráfico
MODEL.DATES.FanchGraph = params.StartDatePlot:params.EndDatePlot;


% DISTRIBUCIÓN
% Percentiles a utilizar (puede variar para agregar más divisiones en la
% distribución: percentiles del 10 al 90 en pasos de 5 por ciento
grilla = params.Grilla;
% del 5 al 95 grilla = 0.05:0.05:0.950001;

% Distribución a utilizar: Normal Estándar (Inversa)
dist = norminv(grilla);
% Cantidad de percentiles
s_dist = size(dist);


%% GRÁFICO
for i = 1:length(params.PlotList)
    
    % Apertura del abanico (simétrica)
    if isempty(params.apertura{i})
       params.apertura{i} = ones(length(MODEL.DATES.Fanchart),1); 
    end
    % Sesgo (asímetría al alza/baja)
    if isempty(params.sesgo{i})
        params.sesgo{i} = ones(1, s_dist(2));
    end

    % Creación de matriz para gráfica y cálculos
    % Error Estándar ajustado
    if length(FMSE.(params.PlotList{i})) == length(params.apertura{i})
        std_adj = FMSE.(params.PlotList{i}).*params.apertura{i};
    else
       error('ERROR: El número de períodos no es consistente entre el FMSE (%d) y del vector de ajuste de sesgo (%d)',...
                      length(FMSE.(params.PlotList{i})), length(params.apertura{i}));
    end
    
    % percentiles ajustados por sesgo
    if length(params.Grilla) == length(params.sesgo{i})
        dist_adj = dist.*params.sesgo{i};
    else
       error('ERROR: El tamaño de la grilla (%d) y del vector de ajuste (%d) de sesgo no es igual',...
              length(params.Grilla), length(params.sesgo{i}));
    end

    % Construcción de matriz
    FAN = zeros(length(MODEL.DATES.FanchGraph),s_dist(2)+1);
    % Primera Columna con serie (historia + pronósticos)
    FAN(:,1) = MODEL.F_pred.(params.PlotList{i}){MODEL.DATES.FanchGraph};
    % Filas 2 a fin de historia, Columnas 2 a última, en períodos de historia (observadas)
    FAN(1:length(MODEL.DATES.FanchGraph(1):MODEL.DATES.hist_end),2:end) = ...
        zeros(length(MODEL.DATES.FanchGraph(1):MODEL.DATES.hist_end),s_dist(2)) + ...
        MODEL.F_pred.(params.PlotList{i}){MODEL.DATES.FanchGraph(1):MODEL.DATES.hist_end}.data;
    % Filas inicio de pronóstico a última, columnas 2 a última, valores de
    % variable en distribución
    FAN(length(MODEL.DATES.FanchGraph(1):MODEL.DATES.hist_end)+1:end,2:end) = ...
        MODEL.F_pred.(params.PlotList{i})(MODEL.DATES.Fanchart)+std_adj.*dist_adj;
    % Cálculo de diferencias para gráfico de stacked areas
    DIFF_FAN = [FAN(:,1), FAN(:,2:end) - FAN(:,1:end-1)];
    DIFF_FAN = tseries(MODEL.DATES.FanchGraph, DIFF_FAN);

    
    %% GRÁFICA

    % Colores
    codigo_color_ref=[1 0 0];
    color_base=0.765*([1 1 1]-codigo_color_ref);
    col_db_fan=size(DIFF_FAN);
    col_db_fan = col_db_fan(2);
    salto=col_db_fan/2-1;
    cambio=color_base./salto;
    codigo_color=codigo_color_ref+color_base;
    inflex=-1;

for ji=1:1:col_db_fan
    if ji~=10
        codigo_color=[codigo_color;codigo_color(ji,:)+inflex*cambio];
    else
        codigo_color=[codigo_color;codigo_color(ji,:)];
        inflex=1;
    end
end
codigo_color=abs(codigo_color);

%% Gráfica
        figure('Position', [1 42.0182 1117.1 776.73]);
        
        h = area(MODEL.DATES.FanchGraph, DIFF_FAN);
        for ii = 1:col_db_fan
            if ii==1
                set(h(ii),'LineStyle','none');
                set(h(ii),'FaceColor','none');
            elseif ii==salto+2
                set(h(ii),'LineWidth',2,'EdgeColor',codigo_color(ii-1,:));
                set(h(ii),'FaceColor',codigo_color(ii-1,:));
            else
                set(h(ii),'LineStyle','none');
                set(h(ii),'FaceColor',codigo_color(ii-1,:));
            end
        end

        hold on
        p = plot(MODEL.DATES.FanchGraph,...
            [MODEL.F_pred.(params.PlotList{i}){MODEL.DATES.FanchGraph},...
             MODEL.F.(params.PlotList{i}){MODEL.DATES.FanchGraph}]);
        set(p(1), 'Color', 'k','Linewidth',2,'LineStyle','--')
        set(p(2), 'Color','k','Linewidth',2)
        
        title(sprintf('%s\nObservado y Pronosticado',MODEL.F_pred.(params.PlotList{i}).comment{1}),...
              'FontSize', 18);
        vline(MODEL.DATES.hist_end, 'LineStyle',':');  
        if strcmp(params.PlotList{i},'D4L_CPI')
            colorarea(MODEL.DATES.FanchGraph,[3 5],[(242/255) (242/255) (242/255)]);
            rng = get(gca,'Xlim');
            a = line([rng(1) rng(2)],[4 4]);set(a,'Color','k','LineStyle','-');
            c = line([rng(1) rng(2)],[3 3]);set(c,'Color','k','LineStyle',':');
            d = line([rng(1) rng(2)],[5 5]);set(d,'Color','k','LineStyle',':'); 
            
            
        else
            SS = get(MODEL.M, 'sstate');
            hline(...
                real(SS.(params.PlotList{i})), ...
                'LineWidth', 0.75, ...
                'LineStyle', ':' ...
                );
        end
        
        hold off
        
                axis on
        
        SimTools.scripts.pausaGuarda(...
            fullfile(params.SavePath, ...
            sprintf("%s.png", params.PlotList{i})), ...
            'AutoSave', params.AutoSave ...
            );


        
        if  isempty(params.Esc)
            MODEL.Esc.v0.Fanchart.(params.PlotList{i}).FAN =  tseries(MODEL.DATES.FanchGraph, FAN);
            MODEL.Esc.v0.Fanchart.(params.PlotList{i}).DIFF_FAN =  DIFF_FAN;
            MODEL.Esc.v0.Fanchart.(params.PlotList{i}).FMSE =  FMSE.(params.PlotList{i});
        else
            MODEL.Esc.(params.Esc).Fanchart.(params.PlotList{i}).FAN =  tseries(MODEL.DATES.FanchGraph, FAN);
            MODEL.Esc.(params.Esc).Fanchart.(params.PlotList{i}).DIFF_FAN =  DIFF_FAN;
            MODEL.Esc.(params.Esc).Fanchart.(params.PlotList{i}).FMSE =  FMSE.(params.PlotList{i});
        end
       
end


end