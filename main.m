%%
%{
QPM v3
%}

clc;
clear all;
PATH.data = genpath('data');
PATH.src = genpath('src');
PATH.temp = genpath('temp');

structfun(@addpath, PATH)

% return
%% Carga de configuraciones generales del corrimiento base (v0)
tic
StartConfig_QPM;

% DEFINIR EL NOMBRE DEL ESCENARIO
scenario_name = 'prueba';
% return

%% Preprocessing
Preprocessing;
disp('Preprocesamiento de Observables: ok');

%% Lectura de Modelo, proceso de filtrado y simulación
MODEL = Sim.sim.read_model(MODEL);

% Filtrado 
[MODEL.MF,MODEL.F] = filter(MODEL.M, MODEL.PreProc.obs,... 
                            MODEL.DATES.hist_start:MODEL.DATES.hist_end, ... 
                            'meanOnly=',true);

% Creación de plan de simulación para escenarios de simulación
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

%% Almacenamiento de Estructura MODEL del mes corriente.
if ~isfolder(fullfile('data','fulldata',MODEL.CORR_DATE))
     mkdir(fullfile('data','fulldata',MODEL.CORR_DATE))
 else
     rmdir(fullfile('data','fulldata',MODEL.CORR_DATE), 's')
     mkdir(fullfile('data','fulldata',MODEL.CORR_DATE))
end 


save(fullfile('data','fulldata',MODEL.CORR_DATE, sprintf("MODEL-%s-QPM.mat", scenario_name)), 'MODEL');
disp('Almacenamiento estructura MODEL: ok');
disp('---- FIN ----');
toc
