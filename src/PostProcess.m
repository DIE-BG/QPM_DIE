%% Variables
pp_list = {'L_MB', 'L_VEL', 'L_CPI_RW', 'L_CPI_RW_Q','L_Z', 'L_GDP', 'L_GDP_RW'};
list_nivel = {'L_S','L_MB'};

% Reconstrucción y pronósticos de PIB de EEUU 
MODEL = rec_GDP_RW(MODEL);
   
% Post-procesamiento
MODEL = PostProcessing(MODEL,...
    'list',pp_list,...
    'list_niv', list_nivel,...
    'Esc',{MODEL.CORR_VER, MODEL.F_pred});

disp('Postprocesamiento: ok');