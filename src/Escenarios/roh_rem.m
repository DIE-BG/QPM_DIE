%{
Se prueba las remesas para distintos valores de Rho para la convergencia
del de la brecha a su estado estacionario. Para este ejercicio se dejó el
estado estacionario en 15% del PIB.
%}

%% Lista de estados rho_rem_gdp a evaluar
rho_rem_gdp = [0.70:0.01:0.99];

%% Generación de modelos a utilizar
% parámetros del QPM
setparam;

% Carga del modelo
for i = 1:length(rho_rem_gdp)
    
    ss = s;
    ss.rho_REM_GDP = rho_rem_gdp(i);
    M = model(MODEL.mod_file_name, 'assign', ss);
    M = sstate(M,'growth=',true,'MaxFunEvals',1000,'display=','off');
    M = solve(M,'error=',true);
    
    fieldName = sprintf('rem_rho_%d', i);
    model_name = sprintf('rem_rho_%d', i);
    
    MODEL.Esc_rem_rho.(fieldName).MODEL = M;
    MODEL.Esc_rem_rho.(fieldName).name = model_name;
    
end


%Filtro de kalman para cada modelo
for i = 1:length(rho_rem_gdp)
    
    fieldName = sprintf('rem_rho_%d', i);
    [MODEL.Esc_rem_rho.(fieldName).MF, MODEL.Esc_rem_rho.(fieldName).F] = filter(MODEL.Esc_rem_rho.(fieldName).MODEL, MODEL.PreProc.obs,... 
                                MODEL.DATES.hist_start:MODEL.DATES.hist_end, ... 
                                'meanOnly=',true);
end

% Pronosticos para cada modelo
fcstrng = MODEL.DATES.pred_start:MODEL.DATES.pred_end;

for i = 1:length(rho_rem_gdp)
    
    fieldName = sprintf('rem_rho_%d', i);
    MODEL.Esc_rem_rho.(fieldName).F_pred = simulate(MODEL.Esc_rem_rho.(fieldName).MF,... Modelo Filtrado
                            MODEL.Esc_rem_rho.(fieldName).F,...Base de datos inicial filtrada
                            fcstrng,... Rango de Simulación
                            'anticipate', false,...
                            'DbOverlay=', true);

end

for i = 1:length(rho_rem_gdp)
    
    fieldName = sprintf('rem_rho_%d', i);
    plot(qq(2024,1):qq(2032,1),...
        MODEL.Esc_rem_rho.(fieldName).F_pred.L_GDP_GAP);
    
    hold on
    
    zeroline();
    
end

MODEL.Esc_rem_rho.rem_rho_15.F_pred.L_GDP_GAP
MODEL.Esc_rem_rho.rem_rho_15.F_pred.rho_REM_GDP
