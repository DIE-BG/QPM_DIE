function plot_diff_shd_dsc(MODEL, varargin)

% plot_diff_shd_dsc realiza las gráficas de descomposición de choques del modelo
% en primeras diferencias.
%{
## Syntax ##

    plot_diff_shd_dsc(MODEL, varargin)

## Input Arguments ##

__`MODEL`__ [ struct ] -
Debe contener al menos la estructura del modelo `MODEL.MF`,
`MODEL.diff_shd_dsc` y `MODEL.DATES`


## Options ##

* SavePath = fullfile(userpath, 'temp') [ `String` ] - Directorio donde
guarda la gráfica.

* Variables = get(MODEL.MF, 'xlist') [ `cell` ] - Nombre de variable a
graficar.

* CloseAll = [ `true`|false ] - Nombre de variable a
graficar.

* 'Esc_add' = [`struct`] - Nombre del escenario, descomposición de shocks y
  y predicciones del escenario.

## Output Arguments ##


## Description ##


## Example ##

%}

% -DIE
% -Octubre 2021

% Parametros opcionales
 p = inputParser;
    addParameter(p, 'SavePath', {});
    addParameter(p, 'Rng', {MODEL.DATES.hist_end-20, MODEL.DATES.hist_end+20});
    addParameter(p, 'Variables', get(MODEL.MF, 'xlist'));
    addParameter(p, 'Esc_add', {}); % libre v0, alterno v1, contrafactual v2...
    addParameter(p, 'CloseAll', true);
    addParameter(p, 'StartDate', {MODEL.DATES.hist_end - 20});
    addParameter(p, 'EndDatePlot', {MODEL.DATES.pred_end});
    addParameter(p, 'MF', {MODEL.MF});
parse(p, varargin{:});
params = p.Results; 

% Verificación y creación del directorio para las gráficas
if ~isfolder(params.SavePath)
    mkdir(params.SavePath)
else
    rmdir(params.SavePath, 's')
    mkdir(params.SavePath)
end

% Variables a descomponer
var_plot = params.Variables;

% Choques a las variables
var_shd = get(params.MF, 'elist');

% Paleta de colores
col = distinguishable_colors(length(var_shd) + 1, ...
    'b', ...
    @(x) colorspace('RGB->Lab',x));

%% calculo de primeras diferencias
% datos
if ~strcmp(params.Esc_add{1}, 'v0')
    full_data_add.shd_dsc = params.Esc_add{2};
    full_data_add.F_pred = params.Esc_add{3};
end

% primeras diferencias
if strcmp(params.Esc_add{1}, 'v0')
    MODEL.diff_shd_dsc = databank.apply(@(x) diff(x),...
                                    MODEL.shd_dsc,...
                                    'InputNames=', params.Variables,...
                                    'RemoveSource=', true);
                                
    MODEL.diff_shd_dsc = MODEL.diff_shd_dsc*params.Variables;
     
else
    full_data_add.diff_shd_dsc = databank.apply(@(x) diff(x),...
                                    full_data_add.shd_dsc,...
                                    'InputNames=', params.Variables,...
                                    'RemoveSource=', true);
                                
    full_data_add.diff_shd_dsc = full_data_add.diff_shd_dsc*params.Variables;
end

%% Graficas
    if strcmp(params.Esc_add, 'v0')

        % El iterador i representa la variable a ser descompuesta
            for i = 1:length(var_plot)

                figure('Position', [1 42.0182 1.6756e+03 825.6000]);
                hold on

                % Barras
                barcon(params.Rng{1}:params.Rng{2},...
                       MODEL.diff_shd_dsc.(var_plot{i}){:, 1:end}, ...
                    'dateFormat=','YYFP', ...
                    'colorMap=',col, ...
                    'evenlySpread=', false); 
                % Líneas
                plot(params.Rng{1}:params.Rng{2},...
                    diff(MODEL.F_pred.(var_plot{i})), ...
                    'w', ...
                    'LineWidth',5);
                plot(params.Rng{1}:params.Rng{2},...
                    diff(MODEL.F_pred.(var_plot{i})), ...
                    'k.-', ...
                    'LineWidth',2, ...
                    'MarkerSize', 20);
                % Línea vertical en el fin de historia
                grfun.vline(MODEL.DATES.hist_end,'timePosition','middle');
                % Título
                set(gca,'FontSize',12);       
                title( ...
                    { ...
                        var_plot{i}; ...
                        "Primera Diferencia" ...
                    }, ...
                    'Interpreter', 'none');
                % Leyendas
                legend(var_shd,'location','northeastoutside','FontSize',11, 'Interpreter', 'none')
                grid on;
                hold off

            % Almacenamiento    
                saveas(gcf, ...
                    fullfile(params.SavePath, ...
                    sprintf("%s_diff_shd_dsc.png", var_plot{i})))
            close
            
            end
            
    else

            % El iterador i representa la variable a ser descompuesta
            for i = 1:length(var_plot)

                figure('Position', [1 42.0182 1.6756e+03 825.6000]);
                hold on

                % Barras
                barcon(params.Rng{1}:params.Rng{2},...
                    full_data_add.diff_shd_dsc.(var_plot{i}){:, 1:end}, ...
                    'dateFormat=','YYFP', ...
                    'colorMap=',col, ...
                    'evenlySpread=', false); 
                % Líneas
                plot(params.Rng{1}:params.Rng{2},...
                    diff(full_data_add.F_pred.(var_plot{i})), ...
                    'w', ...
                    'LineWidth',5);
                plot(params.Rng{1}:params.Rng{2},...
                    diff(full_data_add.F_pred.(var_plot{i})), ...
                    'k.-', ...
                    'LineWidth',2, ...
                    'MarkerSize', 20);
                % Línea vertical en el fin de historia
                grfun.vline(MODEL.DATES.hist_end,'timePosition','middle');
                % Título
                set(gca,'FontSize',12);       
                title( ...
                    { ...
                        var_plot{i}; ...
                        "Primera Diferencia" ...
                    }, ...
                    'Interpreter', 'none');
                % Leyendas
                legend(var_shd,'location','northeastoutside','FontSize',11, 'Interpreter', 'none')
                grid on;
                hold off

            % Almacenamiento    
                saveas(gcf, ...
                    fullfile(params.SavePath, ...
                    sprintf("%s_diff_shd_dsc.png", var_plot{i})))
            close
            end
    end
end
% ----- TEST -----
% params.Variables = get(MODEL.MF, 'xlist');
% i = 1;