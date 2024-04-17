function PostPrLogs(MODEL, varargin)

%% Logaritmos desestacionalizados y tendencia

p = inputParser;
    addParameter(p, 'StartDate', MODEL.DATES.hist_start);
    addParameter(p, 'EndDatePlot', MODEL.DATES.pred_end);
    addParameter(p, 'SavePath', {});
    addParameter(p, 'Esc',{MODEL.CORR_VER, MODEL.F_pred});
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


%%
for rng = 1:length(params.StartDate)
    for var = 1:length(list)
        
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
            'Position', [0, 1 - 0.9, 1, 0.9], ...
            'BackgroundColor', [1, 1, 1] ...
            );
        
        ax = axes(plot_p, 'Units','normalized' ,'Position', [0.1 0.1 0.85 0.8]);
        
        plot(...
            params.StartDate{rng}:params.EndDatePlot{rng}, ...
            MODEL.PostProc.(MODEL.CORR_VER).l_sa.(strcat(list{var},'_sa')),'.-b', ...
            'LineWidth', 2 ...
            );
        
        ax = gca;
        ax.YAxis.Exponent = 0;
        
        hold on
        
        plot(...
            params.StartDate{rng}:params.EndDatePlot{rng}, ...
            MODEL.PostProc.(MODEL.CORR_VER).l_bar.(strcat(list{var},'_bar')),'.-k',...
            'LineWidth', 1.65, ...
            'LineStyle', '--' ...
            );
        
        ax = gca;
        ax.YAxis.Exponent = 0;
        
        hold off
        %Returns handles to the patch and line objects
        chi = get(gca, 'Children');
        %Reverse the stacking order so that the patch overlays the line
        set(gca, 'Children',flipud(chi));
        
        
        
        if ~isempty(params.LegendsNames)
            legend(params.LegendsNames, 'Location', params.LegendLocation);
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
            'Position', [0, 1 - 0.9 - 0.10, 1, 0.10], ...
            'BackgroundColor', [1, 1, 1] ...
            );
        
        data_table = [];
%         if ~isempty(params.FullDataAnt_Name)
            data_table(:, 1) = MODEL.PostProc.(MODEL.CORR_VER).l_bar.(strcat(list{var},'_bar'))(params.TabRange);
            data_table(:, 2) = MODEL.PostProc.(MODEL.CORR_VER).l_sa.(strcat(list{var},'_sa'))(params.TabRange);
            text_Color = [0,0,0 ; 0,0,1];
%         else
%             data_table(:, 1) = MODEL.PostProc.(MODEL.CORR_VER).l_bar.(strcat(list{var},'_bar'))(params.TabRange);
%             text_Color = [1, 1, 1];
%         end
        
        SimTools.scripts.plot_data_table( ...
            params.TabRange, ...
            data_table, ...
            'Parent', table_p, ...
            'SeriesNames', params.LegendsNames, ...
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