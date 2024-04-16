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

%% Post-Procesamiento de variables seleccionadas.
pp_list = {'L_MB', 'L_VEL', 'L_CPI_RW', 'L_CPI_RW_Q'};
list_nivel = {'L_S','L_MB'};

MODEL = PostProcessing(MODEL,...
    'list',pp_list,...
    'list_niv', list_nivel,...
    'Esc',{MODEL.CORR_VER, MODEL.F_pred});
disp('Postprocesamiento: ok');

return
%%

%% Almacenamiento de Estructura MODEL del mes corriente.
save(fullfile('data','fulldata',MODEL.CORR_DATE, sprintf("MODEL-%s.mat", MODEL.CORR_DATE)), 'MODEL');
disp('Almacenamiento estructura MODEL: ok');
disp('---- FIN ----');

