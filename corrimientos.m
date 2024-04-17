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
MODEL.Esc.v0 = MODEL.F_pred;

%% Plots
% Pre-Procesamiento
    % Mensuales
    PreProcPlots_m(MODEL,...
            'Esc_add', {'v0', MODEL_ANT},...
            'tab_range_mm', tab_range_mm);

    % Trimestrales
    PreProcPlots_q(MODEL,...
            'Esc_add', {'v0', MODEL_ANT},...
            'tab_range', tab_range_source_data,...
            'tab_range_mm', tab_range); 
    
% Otras Gráficas
    % Componentes Tipo de cambio real
    tc_real(MODEL,...
        'Esc_add', {'v0', MODEL_ANT.F_pred},...
        'tab_range', tab_range);
    
    % Subplot tcr
    tcr_subplot(MODEL,...
        'Esc_add', {'v0', MODEL_ANT},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.leg_ant, MODEL.leg_act})
    
%% Almacenamiento de Estructura MODEL del mes corriente.
save(fullfile('data','fulldata',MODEL.CORR_DATE, sprintf("MODEL-%s.mat", MODEL.CORR_DATE)), 'MODEL');
disp('Almacenamiento estructura MODEL: ok');
disp('---- FIN ----');

