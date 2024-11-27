%% Modelo y setparam
MODEL.mod_file_name = 'QPM.model';
MODEL.param_file_name = fullfile('src', 'setparam.m');

%% Configuración del corrimiento

MODEL.CORR_VER = 'v0';

MODEL.CORR_DATE     = '2024-11';
MODEL.CORR_DATE_ANT = '2024-09';

% Fechas de fin de historia
MODEL.DATES.hist_end_ant = qq(2024, 2);
MODEL.DATES.hist_end = qq(2024, 3); 
MODEL.DATES.hist_end_mm = mm(2024, 10);

% Nombres de meses para corrimientos automáticos
datestrfn = @(d) datestr(datetime(d, 'InputFormat', 'yyyy-MM'), 'mmmm yyyy', 'local'); 
MODEL.leg_act = datestrfn(MODEL.CORR_DATE);  
MODEL.leg_ant = datestrfn(MODEL.CORR_DATE_ANT); 


%% Otros elementos y fechas
MODEL.data_file_name = fullfile( ...
    'data','corrimientos', MODEL.CORR_DATE, MODEL.CORR_VER, 'data.csv');

MODEL.FULLDATANAME_ACT = fullfile( ...
    'data', 'fulldata', MODEL.CORR_DATE, ...
    sprintf('fulldata_%s_%s.csv', MODEL.CORR_DATE,MODEL.CORR_VER));

MODEL.FULLDATANAME_ANT = fullfile( ...
    'data', 'fulldata', MODEL.CORR_DATE_ANT,...
    sprintf("MODEL-%s.mat", MODEL.CORR_DATE_ANT));

% Configuración de estructura DATES

MODEL.DATES.hist_start = qq(2005, 1);
MODEL.DATES.pred_start = MODEL.DATES.hist_end + 1;
MODEL.DATES.pred_end = MODEL.DATES.hist_end + 30;

% Rango de tablas para gráficos de simulación
tab_range = [MODEL.DATES.hist_end, MODEL.DATES.pred_start:MODEL.DATES.pred_start+5, qq(2026,4)];

% Rango de tablas para gráficos de Pre - procesamiento
% Trimestral
tab_range_source_data = MODEL.DATES.hist_end-8:MODEL.DATES.hist_end;
% Mensual
MODEL.DATES.hist_start_mm = mm(2005,1);
tab_range_mm = MODEL.DATES.hist_end_mm-8:MODEL.DATES.hist_end_mm;

%% Configuración Escenarios
% Nombres de escenarios
MODEL.esc_names = {'Escenario Libre',...v0
                   'Escenario IPEI',...v1
                   'Escenario Tasa Líder',...v2
                   'Escenario Combinado',...v3
                   'Escenario Tasa Externa',...v4
                   'Escenario Anclajes Corto Plazo',...v5
                   'Escenario ss_rem_GDP'};%,...v6


MODEL.Esc.v1.name = MODEL.esc_names{2};
MODEL.Esc.v2.name = MODEL.esc_names{3};
MODEL.Esc.v3.name = MODEL.esc_names{4};
MODEL.Esc.v4.name = MODEL.esc_names{5};
MODEL.Esc.v5.name = MODEL.esc_names{6};
MODEL.Esc.v6.name = MODEL.esc_names{7};


% Colores para escenarios ALTERNOS
MODEL.esc_col = {[0.4660 0.6740 0.1880],...   v1
                 [0.8500 0.3250 0.0980],...    %v2
                 [0.4940 0.1840 0.5560],...    %v3
                 [0.9800 0.0001 0.9500],...    %v4
                 [0.2510 0.1843 0.0156],...    %v5
                 [0.2200 0.9020 1.0000],...    %v6
                 [0.2200 0.9020 1.0000],...    %v7
                 [0.9690 0.5920 0.1490]}; ...  %v8
                      
%% Carga de info mes previo
MODEL_ANT = load(sprintf('MODEL-%s-QPM.mat',MODEL.CORR_DATE_ANT));
MODEL_ANT = MODEL_ANT.MODEL;


%% Lista de variables para post-procesamiento
% Logaritmos desestacionalizados, tendencias y brechas
pp_list = {'L_MB', 'L_VEL', 'L_CPI_RW', 'L_IPEI_Q'};
% Niveles desestacionalizados y tendencias
list_nivel = {'L_S','L_MB'};

% Variables y titulos para gráficas de reconstrucción de nivel
% (PostPrLevels)
list_lev = {'S','MB'};
tit_lev ={{'Tipo de Cambio Nominal (GTQ/USD)'},...
        {'Base Monetaria (Millones de Quetzales)'}};
    
% Lista de gráficas de brechas
list_gaps = {'L_GDP_RW','L_Z', 'L_GDP','L_CPI_RW', 'L_MB', 'L_VEL', 'L_IPEI_Q'};

% Lista de gráficas para descomposición histórica de shocks
sh_list = {'L_GDP_RW_GAP', 'DLA_CPI_RW', 'RS_RW', 'D4L_CPI_NOSUBY','L_GDP_GAP','DLA_CPIXFE', 'DLA_S', 'D4L_MB', 'RS',...
            'D4L_CPI', 'L_Z_GAP', 'D4L_VEL', 'RR', 'D4L_S', 'RMC', 'MCI'};
