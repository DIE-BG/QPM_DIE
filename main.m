%%
%{
QPM
%}

clear all;
PATH.data = genpath('data');
PATH.src = genpath('src');
PATH.temp = genpath('temp');

structfun(@addpath, PATH)

% return
%% Carga de configuraciones generales del corrimiento base (v0)
tic
StartConfig_QPM;

% return

%% Preprocessing
Preprocessing;
disp('Preprocesamiento de Observables: ok');

%% Lectura de Modelo, proceso de filtrado y simulación
MODEL = SimTools.sim.read_model(MODEL);

% Filtrado 
[MODEL.MF,MODEL.F] = filter(MODEL.M, MODEL.PreProc.obs,... 
                            MODEL.DATES.hist_start:MODEL.DATES.hist_end, ... 
                            'meanOnly=',true);

% Creación de plan de simulación con ajustes para EP
PlanSim_v0;
% Simulación
fcstrng = MODEL.DATES.pred_start:MODEL.DATES.pred_end;
MODEL.F_pred = simulate(MODEL.MF,... Modelo Filtrado
                        MODEL.Esc.v0.dbi,...Base de datos inicial filtrada
                        fcstrng,... Rango de Simulación
                        'plan',MODEL.Esc.v0.planSim,... plan de Simulación
                        'anticipate', false,...
                        'DbOverlay=', true);
                    
%% Post-Procesamiento de variables seleccionadas.
PostProcess;

%% Gráficas
do_graphs = true;

if do_graphs == true
    % Pre-processing
    % monthly
    PreProcPlots_m(MODEL,...
        'Esc_add', {'v0', MODEL_ANT},...
        'tab_range_mm', tab_range_mm);
    
    % quarterly
    PreProcPlots_q(MODEL,...
        'Esc_add', {'v0', MODEL_ANT},...
        'tab_range', tab_range_source_data)
    % Graficas de simulación
    simPlots(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, MODEL.CORR_VER,'prediction_compared'),...}
        'Esc_add', {MODEL.CORR_DATE_ANT, MODEL_ANT.F_pred,[]},...
        'PlotList', get(MODEL.MF, 'xlist'),...
        'LegendsNames',{MODEL.leg_ant, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Post Procesamiento
    % Logaritmo/tendencia vs corr. Anterior
    PostPrLogsComp(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE,'v0', 'PostProcessing'),...
        'Esc',{'v0', MODEL.PostProc},...
        'Esc_add', {'v0', MODEL_ANT.PostProc.v0,[]},...
        'PlotList', pp_list,...
        'LegendsNames',{MODEL.leg_ant, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Variables en Niveles originales
    PostPrLevels(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, MODEL.CORR_VER,'PostProcessing'),...
        'Esc_add', {MODEL.CORR_DATE_ANT, MODEL_ANT.PostProc.v0, []},...
        'PlotList', list_lev,...
        'Titles', tit_lev,...
        'LegendsNames',{MODEL.leg_ant, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Brechas
    PostPrGaps(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, MODEL.CORR_VER,'PostProcessing'),...
        'Esc_add', {MODEL.CORR_DATE_ANT, MODEL_ANT.PostProc.v0,[]},...
        'PlotList', list_gaps,...
        'LegendsNames',{MODEL.leg_ant, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Descomposición de choques para variables seleccionadas
    Desc_shocks;
    
    % Descomposición de chosques para variables seleccionadas
    % (Primeras diferencias)
    plot_diff_shd_dsc(MODEL,...
                  'variables', sh_list,...
                  'SavePath', fullfile('plots', MODEL.CORR_DATE, MODEL.CORR_VER, 'Shock_dec', 'diff'),...
                  'Esc_add', {'v0'});
    
    % Contribuciones
    contributions(MODEL,...
                  'Esc_add', {'v0', MODEL_ANT});
    diff_contributions(MODEL,...
                  'Esc_add', {'v0', MODEL_ANT});
              
    % Real exchange rate (subplot)
    tcr_subplot(MODEL,...
        'Esc_add', {'v0', MODEL_ANT.F_pred, []},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.leg_ant, MODEL.leg_act});
    
    % Money Velocity (subplot)
    vel_subplot(MODEL,...
        'tab_range', tab_range,...
        'Esc_add', {'v0', MODEL_ANT.PostProc, []},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.leg_ant, MODEL.leg_act});
    
     MODEL = fanchart(MODEL, MODEL_ANT,...
                    'EndDatePlot',qq(2025,4),...
                    'apertura', {[0.25 0.25 0.25 0.25 0.25 0.25 0.25]'},... Vector Columna
                    'sesgo', {[1 1 1 1,... Percentiles por debajo de la mediana
                               1,... percentil 50
                               1.5 1.5 1.5 1.5],... percentiles por arriba de la mediana
                               }...Vector fila
                    );  
    
end

%% fulldata libre

if ~isfolder(fullfile('data', 'fulldata', MODEL.CORR_DATE))
    mkdir(fullfile('data', 'fulldata', MODEL.CORR_DATE))
end

%Creación de fulldata para escenario IPEI del SVAR-50-4B
databank.toCSV(MODEL.F_pred*{'L_CPI_RW','DLA_CPI_RW', 'D4L_CPI_RW' }, fullfile('data', 'fulldata', MODEL.CORR_DATE, 'fulldata_QPM.csv'), Inf);

%Creación de fulldata completo del escenario libre QPM
databank.toCSV(MODEL.F_pred, fullfile('data', 'fulldata', MODEL.CORR_DATE, 'fulldata.csv'), Inf);

% return; %Correr hasta esta linea la primera vez para Generar los anclajes necesarios para el Escneario IPEI del SVAR
%% Escenarios alternos
esc_alt = true;
graph_esc = true;

if esc_alt == true
    v1_IPEI;
    v2_CP1;
    v3_Comb;
    v4_Istar;
    v5_Anclaje_mm;
end


%% Presentación
prs = true;
if prs == true
    presentacion;
end

%% Almacenamiento de Estructura MODEL del mes corriente.
save(fullfile('data','fulldata',MODEL.CORR_DATE, sprintf("MODEL-%s-QPM.mat", MODEL.CORR_DATE)), 'MODEL');
disp('Almacenamiento estructura MODEL: ok');
disp('---- FIN ----');
toc
