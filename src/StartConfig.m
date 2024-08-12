%% Modelo y setparam
MODEL.mod_file_name = 'QPM.model';
MODEL.param_file_name = 'setparam.m';

%% Configuración del corrimiento

MODEL.CORR_VER = 'v0';

MODEL.CORR_DATE = '2024-05';
MODEL.CORR_DATE_ANT = '2024-04';

MODEL.leg_act = 'Mayo 2024';  
MODEL.leg_ant = 'Abril 2024'; 

% Fechas de fin de historia
MODEL.DATES.hist_end_ant = qq(2024, 1);
MODEL.DATES.hist_end = qq(2024, 1);%qq(2023,4);%
MODEL.DATES.hist_end_mm = mm(2024, 04);

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
tab_range = [MODEL.DATES.hist_end, MODEL.DATES.pred_start:MODEL.DATES.pred_start+5, qq(2025,4)];
arrow_range = qq(2023,4):4:qq(2025,4);
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
                   'Escenario Combinado'};%,...v3
% Colores para escenarios ALTERNOS
MODEL.esc_col = {[0.4660 0.6740 0.1880],...   v1
                 [0.8500 0.3250 0.0980],...    %v2           
                 [0.4940 0.1840 0.5560]}; %v3
                      
%% Carga de info mes previo
MODEL_ANT = load(fullfile('data','fulldata',MODEL.CORR_DATE_ANT,sprintf('MODEL-%s.mat',MODEL.CORR_DATE_ANT)));
MODEL_ANT = MODEL_ANT.MODEL;


%% Lista de variables para post-procesamiento
% Logaritmos desestacionalizados, tendencias y brechas
pp_list = {'L_MB', 'L_VEL', 'L_CPI_RW', 'L_CPI_RW_Q'};
% Niveles desestacinoalizados y tendencias
list_nivel = {'L_S','L_MB'};

% Variables y titulos para gráficas de reconstrucción de nivel
% (PostPrLevels)
list_lev = {'S','MB'};
tit_lev ={{'Tipo de Cambio Nominal (GTQ/USD)'},...
        {'Base Monetaria (Millones de Quetzales)'}};
    
% Lista de gráficas de brechas
list_gaps = {'L_GDP_RW','L_Z', 'L_GDP','L_CPI_RW', 'L_MB', 'L_VEL', 'L_CPI_RW_Q'};
