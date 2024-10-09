function PostPrGaps(MODEL, varargin)

%{
% Brechas vrs mes anterior u otro escenario.
    Genera las gráficas de los niveles desestacionalizados. Pueden compararse,
    o no, con otro escenario o con otro corrimiento.
{
## Syntax ##

    MODEL = simPlots(MODEL, varargin)

## Input Arguments ##

__`MODEL`__ [ struct ] -
Debe contener al menos la estructura con los resultados del proceso de
simulación MODEL.F_pred.

* 'StartDate' = {} [ `Cell` ] - fechas de inicio del plot (pueden ser una o mas).

* 'EndDatePlot' = {} [ `Cell` ] - fechas de fin del plot (pueden ser una o mas).

* 'Esc_add' = {}  [ `Cell` ] - Escenario adicional a plotear. Cell array 
    con dos elementos: (1) Versión o fecha del escenario adicional y (2)
    Base de datos con los pronósticos del modelo (ambas estructuras deben
    tener los mismos campos por compatibilidad.

* 'PlotList' = {} [ `Cell` ] - Lista de brechas a plotear.

- DIE
- Marzo 2024
- MJGM/JGOR
%}

p = inputParser;
    addParameter(p, 'StartDate', MODEL.DATES.hist_start);
    addParameter(p, 'EndDatePlot', MODEL.DATES.pred_end);
    addParameter(p, 'SavePath', fullfile(userpath, 'temp'));
    addParameter(p, 'Esc_add', {});
    addParameter(p, 'PlotList', {});
    addParameter(p, 'LegendsNames', {});
    addParameter(p, 'LegendLocation', 'best');
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
%%

if isempty(params.SavePath)
    params.SavePath = fullfile('plots', MODEL.CORR_DATE, params.Esc{1}, 'PostProcessing');
end

if ~isfolder(params.SavePath)
    mkdir(params.SavePath)
end  

%% Carga de base de datos adicional


% strcmp(params.Esc_add{1}, MODEL.CORR_DATE_ANT)
if ~isempty(params.Esc_add)
    PostProcAdd = params.Esc_add{2}; 
    esc_col = params.Esc_add{3};
end
if isempty(params.Esc_add{3})
    esc_col = 	[1 0 0];
end

%%

list = params.PlotList;

for rng = 1:length(params.StartDate)
    
    for var = 1:length(list)
        
        
        figure;
        
        set(gcf, ...
            'defaultaxesfontsize',12, ...
            'Position', [1 42.0182 1117.1 776.73] ...
            );
        main_p = uipanel('Units','normalized');
        
        % ----- Panel de gráfica -----
        plot_p = uipanel( ...
            main_p, ...
            'Position', [0, 1 - 0.9, 1, 0.9], ...
            'BackgroundColor', [1, 1, 1] ...
            );
        
        ax = axes(plot_p, 'Units','normalized' ,'Position', [0.1 0.1 0.85 0.8]);
        
        plot(...
            params.StartDate{rng}:params.EndDatePlot{rng}, ...
            MODEL.PostProc.v0.l_gap.(strcat(list{var},'_GAP')),'.-b', ...
            'MarkerSize', 17, ...
            'LineWidth', 2 ...
            );
        ax = gca;
        ax.YAxis.Exponent = 0;
        hold on
        
        plot(...
            params.StartDate{rng}:params.EndDatePlot{rng}, ...
            PostProcAdd.l_gap.(strcat(list{var},'_GAP')),...
            'Color',esc_col,...
            'LineWidth', 1.65, ...
            'LineStyle', '--' ...
            );
        
        ax = gca;
        ax.YAxis.Exponent = 0;
        zeroline;
        hold off
        %Returns handles to the patch and line objects
        chi = get(gca, 'Children');
        %Reverse the stacking order so that the patch overlays the line
        set(gca, 'Children',flipud(chi));
        
        
        
        if ~isempty(params.LegendsNames)
            legend(params.LegendsNames, 'Location', params.LegendLocation,...
                   'Interpreter','none');
        end
        
        if strcmp(list{var}, 'L_GDP_RW')
            title( ...
            sprintf(...
            '%s \n(Brecha) \n%s - %s', ...
            'Producto Externo (EEUU)', ...
            dat2char(params.StartDate{rng}), ...
            dat2char(params.EndDatePlot{rng})...
            ) ,...
            'Interpreter','none'...
            );
        else
        title( ...
            sprintf(...
            '%s \n(Brecha) \n%s - %s', ...
            MODEL.F_pred.(list{var}).Comment{:}(1:end-11), ...
            dat2char(params.StartDate{rng}), ...
            dat2char(params.EndDatePlot{rng})...
            ) ,...
            'Interpreter','none'...
            );
        end
        highlight(params.StartDate{rng}:MODEL.DATES.hist_end);
        
        vline( ...
            MODEL.DATES.hist_end,...
            'LineWidth', 1, ...
            'LineStyle', '-.' ...
            );
        
        
        % ----- Panel de Tabla -----
        table_p = uipanel( ...
            main_p, ...
            'Position', [0, 1 - 0.9 - 0.10, 1, 0.10], ...
            'BackgroundColor', [1, 1, 1] ...
            );
        
        data_table = [];
            data_table(:, 1) = PostProcAdd.l_gap.(strcat(list{var},'_GAP'))(params.TabRange);
            data_table(:, 2) = MODEL.PostProc.v0.l_gap.(strcat(list{var},'_GAP'))(params.TabRange);
            text_Color = [esc_col ; 0,0,1];
        
        Sim.scripts.plot_data_table( ...
            params.TabRange, ...
            data_table, ...
            'Parent', table_p, ...
            'SeriesNames', params.LegendsNames, ...
            'TextColor', text_Color, ...
            'FontSize', 9 ...
            )
        axis on
        
        if rng == 2
            save_name = sprintf("%s_GAP_short.png", list{var});
        else
            save_name = sprintf("%s_GAP.png", list{var});
        end
        %
        Sim.scripts.pausaGuarda(...
            fullfile(params.SavePath, ...
            save_name), ...
            'AutoSave', params.AutoSave ...
            );
    end
    if params.CloseAll
        % close all
    end
end


end