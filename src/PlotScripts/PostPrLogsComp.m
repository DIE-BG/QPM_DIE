function PostPrLogsComp(MODEL, varargin)

%{
% Logaritmos desestacionalizados y tendencia
    Genera las gráficas de los logaritmos desestacionalizados y su tendencia. 
    Pueden compararse, o no, con otro escenario o con otro corrimiento.
{
## Syntax ##

    MODEL = simPlots(MODEL, varargin)

## Input Arguments ##

__`MODEL`__ [ struct ] -
Debe contener al menos la estructura con los resultados del proceso de
simulación MODEL.F_pred.

* 'StartDate' = {} [ `Cell` ] - fechas de inicio del plot (pueden ser una o mas).

* 'EndDatePlot' = {} [ `Cell` ] - fechas de fin del plot (pueden ser una o mas).

* 'Esc' = {}  [ `Cell` ] - Escenario principal a plotear. 
    (usualmente escenario libre del corrimiento en curso). Cell array con  
    dos elementos: (1) Versión del escenario principal (v0) y (2)
    Base de datos con los pronósticos del modelo (ambas estructuras deben
    tener los mismos campos por compatibilidad.

* 'Esc_add' = {}  [ `Cell` ] - Escenario adicional a plotear. Cell array 
    con dos elementos: (1) Versión o fecha del escenario adicional y (2)
    Base de datos con los pronósticos del modelo (ambas estructuras deben
    tener los mismos campos por compatibilidad.

* 'PlotList' = {} [ `Cell` ] - Lista de variables a plotear. Compatible con
    lista de logaritmos que entra en función PostProcessing

- DIE
- Marzo 2024
- MJGM/JGOR
%}
p = inputParser;
    addParameter(p, 'StartDate', MODEL.DATES.hist_start);
    addParameter(p, 'EndDatePlot', MODEL.DATES.pred_end);
    addParameter(p, 'SavePath', {});
    addParameter(p, 'Esc',{MODEL.CORR_VER, MODEL.PostProc.(MODEL.CORR_VER)});
    addParameter(p, 'Esc_add', {});
    addParameter(p, 'PlotList', {});
    addParameter(p, 'LegendsNames', {});
    addParameter(p, 'LegendLocation', 'Best');
    addParameter(p, 'PlotSSLine', true);
    addParameter(p, 'PlotAnnotations', true);
    addParameter(p, 'AnnotationXAdjustment', 0);
    addParameter(p, 'AnnotationYAdjustment', 0);
    addParameter(p, 'AnnoRange', qq(2022,4):4:qq(2024,4));
    addParameter(p, 'TabRange', qq(2021,4):4:qq(2024,4));   
    addParameter(p, 'CloseAll', true);    
    addParameter(p, 'AutoSave', true);
 
parse(p, varargin{:});
params = p.Results;

list = params.PlotList;

%% Limpieza y creación de folders
%SavePath
if isempty(params.SavePath)
    params.SavePath = fullfile('plots', MODEL.CORR_DATE, params.Esc{1}, 'PostProcessing');
end

if ~isfolder(params.SavePath)
    mkdir(params.SavePath)
else
    rmdir(params.SavePath, 's')
    mkdir(params.SavePath)
end  

if ~isempty(params.Esc_add)
    full_data_add = params.Esc_add{2};
end

%%
for rng = 1:length(params.StartDate)
    for var = 1:length(list)
    if ~isempty(params.EndDatePlot{rng})
        full_data_add_temp = dbclip(full_data_add, params.StartDate{rng}:params.EndDatePlot{rng});
%         F_pred_temp = dbclip(MODEL.F_pred, params.StartDate{rng}:params.EndDatePlot{rng});
    end
        
        %%
        figure;
        
        set(gcf, ...
            'defaultaxesfontsize',12, ...
            'Position', [1 42.0182 1117.1 776.73] ...
            );
        main_p = uipanel('Units','normalized');
        
        % ----- Panel de gráfica -----
        plot_p = uipanel( ...
            main_p, ...
            'Position', [0, 1 - 0.8, 1, 0.8], ...
            'BackgroundColor', [1, 1, 1] ...
            );
        
        ax = axes(plot_p, 'Units','normalized' ,'Position', [0.1 0.1 0.85 0.8]);
        
        plot(...
            params.StartDate{rng}:params.EndDatePlot{rng}, ...
            MODEL.PostProc.(MODEL.CORR_VER).l_sa.(strcat(list{var},'_SA')),'.-b', ...
            'LineWidth', 2, ...
            'DisplayName',sprintf('Logaritmo - %s', params.LegendsNames{2})...
            );
        
        ax = gca;
        ax.YAxis.Exponent = 0;
        
        hold on
        
        plot(...
            params.StartDate{rng}:params.EndDatePlot{rng}, ...
            MODEL.PostProc.(MODEL.CORR_VER).l_bar.(strcat(list{var},'_BAR')),'.-b',...
            'LineWidth', 1.65, ...
            'LineStyle', '--', ...
            'DisplayName',sprintf('Tendencia - %s', params.LegendsNames{2})...
            );
        
        
        
        if ~isempty(params.Esc_add)
        plot(...
            params.StartDate{rng}:params.EndDatePlot{rng}, ...
            full_data_add_temp.l_sa.(strcat(list{var},'_SA')),'.-r', ...
            'LineWidth', 2, ...
            'DisplayName',sprintf('Logaritmo - %s', params.LegendsNames{1})...
            );
        plot(...
            params.StartDate{rng}:params.EndDatePlot{rng}, ...
            full_data_add_temp.l_bar.(strcat(list{var},'_BAR')),'.-r',...
            'LineWidth', 1.65, ...
            'LineStyle', '--', ...
            'DisplayName',sprintf('Tendencia - %s', params.LegendsNames{1})...
            );        
        
            
        end
        
        ax = gca;
        ax.YAxis.Exponent = 0;
        hold off
        %Returns handles to the patch and line objects
        chi = get(gca, 'Children');
        %Reverse the stacking order so that the patch overlays the line
        set(gca, 'Children',flipud(chi));
        
        
        
        if ~isempty(params.LegendsNames)
            legend('Location', params.LegendLocation);
        end
        
        
        if isempty(MODEL.F_pred.(list{var}).comment)
            temp_title = {};
        else
            temp_title = MODEL.F_pred.(list{var}).comment{1};
        end
        
        title( ...
            sprintf(...
            '%s \n %s - %s', ...
            temp_title, ...
            dat2char(params.StartDate{rng}), ...
            dat2char(params.EndDatePlot{rng})...
            ) ,...
            'Interpreter','none'...
            );
        
        highlight(params.StartDate{rng}:MODEL.DATES.hist_end);
        
        vline( ...
            MODEL.DATES.hist_end,...
            'LineWidth', 1, ...
            'LineStyle', '-.' ...
            );
        
        
        % ----- Panel de Tabla -----
        table_p = uipanel( ...
            main_p, ...
            'Position', [0, 1 - 0.8 - 0.20, 1, 0.20], ...
            'BackgroundColor', [1, 1, 1] ...
            );
        
        data_table = [];
        if isempty(params.Esc_add)
            data_table(:, 1) = MODEL.PostProc.(MODEL.CORR_VER).l_bar.(strcat(list{var},'_BAR'))(params.TabRange);
            data_table(:, 2) = MODEL.PostProc.(MODEL.CORR_VER).l_sa.(strcat(list{var},'_SA'))(params.TabRange);
            text_Color = [0,0,0 ; 0,0,1];
        else
            data_table(:, 1) = full_data_add_temp.l_bar.(strcat(list{var},'_BAR'))(params.TabRange);
            data_table(:, 2) = full_data_add_temp.l_sa.(strcat(list{var},'_SA'))(params.TabRange);
            data_table(:, 3) = MODEL.PostProc.(MODEL.CORR_VER).l_bar.(strcat(list{var},'_BAR'))(params.TabRange);
            data_table(:, 4) = MODEL.PostProc.(MODEL.CORR_VER).l_sa.(strcat(list{var},'_SA'))(params.TabRange);            
            text_Color = [1,0,0;  1,0,0; 0,0,1 ;0,0,1];
        end
        
        SimTools.scripts.plot_data_table( ...
            params.TabRange, ...
            data_table, ...
            'Parent', table_p, ...
            'SeriesNames', {sprintf('Tendencia - %s', params.LegendsNames{1}),...
                            sprintf('Logaritmo - %s', params.LegendsNames{1}),...
                            sprintf('Tendencia - %s', params.LegendsNames{2}),...
                            sprintf('Logaritmo - %s', params.LegendsNames{2})},...
            'TextColor', text_Color, ...
            'FontSize', 9 ...
            )
        axis on
        
        % if var < 10
        %     num = sprintf("0%i", var);
        % else
        %     num = sprintf('%i', var);
        % end
        
        
        if rng == 2
            save_name = sprintf("%s_short.png", list{var});
        else
            save_name = sprintf("%s.png", list{var});
        end
        
        SimTools.scripts.pausaGuarda(...
            fullfile(params.SavePath, ...
            save_name), ...
            'AutoSave', params.AutoSave ...
            );
    end
    if params.CloseAll
        close all
    end
end

end