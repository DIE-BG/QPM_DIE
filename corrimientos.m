%%
%{
QPM
%}

% clear all;
PATH.data = genpath('data');
PATH.src = genpath('src');
PATH.temp = genpath('temp');

structfun(@addpath, PATH)

%% Carga de configuraciones generales del corrimiento base (v0)
tic
StartConfig;

%% Preprocessing
Preprocessing;
disp('Preprocesamiento de Observables: ok');

%% Lectura de Modelo, proceso de filtrado y simulación
MODEL = SimTools.sim.read_model(MODEL);
% Filtrado 
[MODEL.MF,MODEL.F] = filter(MODEL.M, MODEL.PreProc.obs,... 
                            MODEL.DATES.hist_start:MODEL.DATES.hist_end, ... 
                            'meanOnly=',true);
fcstrng = MODEL.DATES.pred_start:MODEL.DATES.pred_end;
MODEL.F_pred = simulate(MODEL.MF, MODEL.F, fcstrng, 'anticipate', false, 'DbOverlay=', true);

%% Post-Procesamiento de variables seleccionadas.
pp_list = {'L_MB', 'L_VEL', 'L_CPI_RW', 'L_CPI_RW_Q'};
list_nivel = {'L_S','L_MB'};

MODEL = PostProcessing(MODEL,...
    'list',pp_list,...
    'list_niv', list_nivel,...
    'Esc',{MODEL.CORR_VER, MODEL.F_pred});
disp('Postprocesamiento: ok');

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
        'Esc_add', {MODEL.CORR_DATE_ANT, MODEL_ANT.F_pred},...
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
        'Esc_add', {'v0', MODEL_ANT.PostProc.v0},...
        'PlotList', pp_list,...
        'LegendsNames',{MODEL.leg_ant, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Variables en Niveles originales
    PostPrLevels(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, MODEL.CORR_VER,'PostProcessing'),...
        'Esc_add', {MODEL.CORR_DATE_ANT, MODEL_ANT.PostProc.v0},...
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
        'Esc_add', {MODEL.CORR_DATE_ANT, MODEL_ANT.PostProc.v0},...
        'PlotList', list_gaps,...
        'LegendsNames',{MODEL.leg_ant, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Descomposición de choques para variables seleccionadas
    Desc_shocks

    
    % Real exchange rate
    tc_real(MODEL,...
        'Esc_add', {'v0', MODEL_ANT.F_pred},...
        'tab_range', tab_range);
    
    % Real exchange rate (subplot)
    tcr_subplot(MODEL,...
        'Esc_add', {'v0', MODEL_ANT},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.leg_ant, MODEL.leg_act});
    
    % Money Velocity (subplot)
    vel_subplot(MODEL,...
        'tab_range', tab_range,...
        'Esc_add', {'v0', MODEL_ANT.PostProc},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.leg_ant, MODEL.leg_act});
     
    % shock decomposition
    Desc_shocks;
    
end

%%

%% Presentación
prs = true;
if prs == true
    presentacion;
end

%% Almacenamiento de Estructura MODEL del mes corriente.
save(fullfile('data','fulldata',MODEL.CORR_DATE, sprintf("MODEL-%s.mat", MODEL.CORR_DATE)), 'MODEL');
disp('Almacenamiento estructura MODEL: ok');
disp('---- FIN ----');
toc
