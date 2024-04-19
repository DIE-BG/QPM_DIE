function PreProcPlots_m(MODEL, varargin)
%{
Graficas con los datos fuente (Frecuencia mensual)
    - Producto externo
    - Producto doméstico
    - Tasa de interes de fondos federales
    - Indices y tasas de variacion de precios de importaciones y
    exportaciones
    - Evolucion del ponderador alpha
    - Indices de Inflación (total y subyacente)
    - Tipo de cambio
    - Remesas
    - Base Monetaria

__`MODEL`__ [ struct ] -
Debe contener al menos la estructura con los resultados del 
pre-procesamiento de los datos.

* 'StartDate' = {} [ `Cell` ] - Fechas de inicio del plot (Puede ser una o mas).
* 'EndDatePlot' = {} [ `Cell` ] - fechas de fin del plot (pueden ser una o mas).

* 'Esc_add' = {}  [ `Cell` ] - Escenario adicional a plotear. Cell array 
    con dos elementos: (1) Versión o fecha del escenario adicional y (2)
    Base de datos con los datos historicos del corrimiento anterior.

Departamento de Investigaciones Económicas - 2024.
MJGM/JGOR
%}

p = inputParser;
addParameter(p, 'StartDate', {MODEL.DATES.hist_start_mm, MODEL.DATES.hist_end_mm - 60});
addParameter(p, 'EndDatePlot', {MODEL.DATES.hist_end_mm});
addParameter(p, 'SavePath', {});
addParameter(p, 'Esc_add', {}); % libre v0, alterno v1, contrafactual v2...
addParameter(p, 'tab_range', {});
addParameter(p, 'tab_range_mm', {});
addParameter(p, 'TabRange', qq(2021,4):4:qq(2024,4));

parse(p, varargin{:});
params = p.Results;

%% Limpieza y creación de folders
% Verificación y creación del directorio para las gráficas
if isempty(params.SavePath)
    params.SavePath = fullfile('plots', MODEL.CORR_DATE, params.Esc_add{1}, 'PreProcessing');
end

if ~isfolder(params.SavePath)
    mkdir(params.SavePath)
else
    rmdir(params.SavePath, 's')
    mkdir(params.SavePath)
end

%% Carga de base de datos mes anterior
full_data_add = params.Esc_add{2};

%% Leyendas
params.LegendsNames = {MODEL.leg_act, MODEL.leg_ant};

%% Evolución alpha
for rng = params.StartDate

    set( ...
        gcf, ...
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
    
    
    plt = plot(rng{1}:params.EndDatePlot{1},...
        [MODEL.PreProc.monthly.A, MODEL.PreProc.monthly.A_prom], ...
        'LineWidth', 1.25);
    
    % Colores
    % mes con mes
    set(plt(1), 'color', [0 0 0], 'LineWidth', 1.25);
    % promedio movil 12 meses
    set(plt(2), 'color', [0 0.498 0],'LineWidth', 1.25, 'Marker', '.', 'MarkerSize', 12);
    
    % leyenda
    legend({'Evolución mensual', 'Promedio Movil (12 meses)'}, 'Location','best');
    
    %Titulos
    title('Importancia de EEUU en importaciones de GT',...
        strcat(dat2char(rng{1}),...
        '-',dat2char(params.EndDatePlot{1})));
    
    % panel
    table_p = uipanel( ...
        main_p, ...
        'Position', [0, 1 - 0.8 - 0.20, 1, 0.20], ...
        'BackgroundColor', [1, 1, 1] ...
        );
    
    data_table = [];
    data_table(:, 1) = MODEL.PreProc.monthly.A(params.tab_range_mm);
    data_table(:, 2) = MODEL.PreProc.monthly.A_prom(params.tab_range_mm);
    text_Color = [[0 0 0]; [0 0.498 0]];
    
    
    SimTools.scripts.plot_data_table( ...
        params.tab_range_mm, ...
        data_table, ...
        'Parent', table_p, ...
        'SeriesNames', {'Evolución mensual', 'Promedio Movil 12 meses'}, ...
        'TextColor', text_Color, ...
        'ColNameWidth', 0.23, ...
        'FontSize', 11 ...
        );
    
    axis on;
    
    % Guardamos la gráfica
    if rng{1} == params.StartDate{1}
        SimTools.scripts.pausaGuarda(...
            fullfile(params.SavePath,...
            'alpha.png'), ...
            'AutoSave', true);
        
    elseif rng{1} == params.StartDate{2}
        SimTools.scripts.pausaGuarda(...
            fullfile(params.SavePath,...
            'alpha_short.png'), ...
            'AutoSave', true);
    end
    
end

%% Resto de variables en frecuencia mensual
names = dbnames(MODEL.PreProc.monthly);

for rng = params.StartDate
    for i = 1:length(names)
        set( ...
            gcf, ...
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
        
        % Corrimiento actual
        plot(rng{1}:params.EndDatePlot{1},...
            full_data_add.PreProc.monthly.(names{i}),'.-r',...
            'MarkerSize', 12, ...
            'LineWidth', 1.5 ...
            );

        ax = gca;
        ax.YAxis.Exponent = 0;
        
        hold on
        
        % Corrimiento anterior
        plot(rng{1}:params.EndDatePlot{1},...
            MODEL.PreProc.monthly.(names{i}),'.-b',...
            'MarkerSize', 12, ...
            'LineWidth', 1.5 ...
            );
        
        ax = gca;
        ax.YAxis.Exponent = 0;
        
        hold off
        
        if startsWith(names{i}, ["DLA_", "D4L_"]) || endsWith(names{i}, ["mom", "yoy"])
            zeroline()
        end
        
        %Returns handles to the patch and line objects
        chi = get(gca, 'Children');
        %Reverse the stacking order so that the patch overlays the line
        set(gca, 'Children',flipud(chi));
        
        % leyenda
        legend(params.LegendsNames, 'Location','best');

        %Titulos
        if isempty(strmatch(names{i}, {'CPI_mom', 'CPI_yoy',...
                                       'CPIXFE_mom', 'CPIXFE_yoy',...
                                       'CPINOSUBY_mom', 'CPINOSUBY_yoy'},...
                                       'exact'))
        
        title(MODEL.PreProc.monthly.(names{i}).Comment{1},...
            strcat(dat2char(rng{1}),...
            '-',dat2char(params.EndDatePlot{1})));
        
        else
                       
        title({MODEL.PreProc.monthly.(names{i}).Comment{1},...
            MODEL.PreProc.monthly.(names{i}).Caption},...
            strcat(dat2char(rng{1}),...
                '-',dat2char(params.EndDatePlot{1})),...
                'Fontsize', 11);
            
        end

        % panel
        table_p = uipanel( ...
            main_p, ...
            'Position', [0, 1 - 0.8 - 0.20, 1, 0.20], ...
            'BackgroundColor', [1, 1, 1] ...
            );

        data_table = [];
        data_table(:, 1) = MODEL.PreProc.monthly.(names{i})(params.tab_range_mm);
        data_table(:, 2) = full_data_add.PreProc.monthly.(names{i})(params.tab_range_mm);
        text_Color = [[0 0 1]; [1 0 0]];


        SimTools.scripts.plot_data_table( ...
            params.tab_range_mm, ...
            data_table, ...
            'Parent', table_p, ...
            'SeriesNames', {params.LegendsNames{1}, params.LegendsNames{2}}, ...
            'TextColor', text_Color, ...
            'ColNameWidth', 0.23, ...
            'FontSize', 10 ...
            );

        axis on;

        % Guardamos la gráfica
        if rng{1} == params.StartDate{1}
            SimTools.scripts.pausaGuarda(...
                fullfile(params.SavePath,...
                sprintf("m_%s.png", names{i})), ...
                'AutoSave', true);

        elseif rng{1} == params.StartDate{2}
            SimTools.scripts.pausaGuarda(...
                fullfile(params.SavePath,...
                sprintf("m_%s_short.png", names{i})), ...
                'AutoSave', true);
        end
    end
end




end
