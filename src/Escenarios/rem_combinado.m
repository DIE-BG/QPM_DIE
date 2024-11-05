%{
Prueba combinando un estado estacionario de remesas más alto y un
rho_rem_gdp mas bajo
%}

% parámetros del QPM
setparam;
s.ss_REM_GDP = 18.8;
s.rho_REM_GDP = 0.90;

% Modelo
M = model(MODEL.mod_file_name, 'assign', s);
M = sstate(M,'growth=',true,'MaxFunEvals',1000,'display=','off');
M = solve(M,'error=',true);
    
MODEL.MODEL_REM = M;

% Filtro de kalman
[MODEL.MF_REM, MODEL.F_REM] = filter(MODEL.MODEL_REM, MODEL.PreProc.obs,... 
                                MODEL.DATES.hist_start:MODEL.DATES.hist_end, ... 
                                'meanOnly=',true);

% Pronosticos
% Simulación
fcstrng = MODEL.DATES.pred_start:MODEL.DATES.pred_end;
MODEL.F_REM_pred = simulate(MODEL.MF_REM,... Modelo Filtrado
                        MODEL.F_REM,...Base de datos inicial filtrada
                        fcstrng,... Rango de Simulación
                        'anticipate', false,...
                        'DbOverlay=', true);

