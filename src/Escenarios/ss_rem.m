%{
Creación de escenarios modificando el estados estacionario de las remesas
de 15 a 19 en pasos de 0.10 para encontrar una brecha que cierre después de
dos años pero no mayor a 5
%}

%% Lista de estados estacionarios a evaluar
ss_rem_gdp = [17:0.1:19];

%% Generación de modelos a utilizar
% parámetros del QPM
setparam;

% Carga del modelo
for i = 1:length(ss_rem_gdp)
    
    ss = s;
    ss.ss_REM_GDP = ss_rem_gdp(i);
    M = model(MODEL.mod_file_name, 'assign', ss);
    M = sstate(M,'growth=',true,'MaxFunEvals',1000,'display=','off');
    M = solve(M,'error=',true);
    
    fieldName = sprintf('rem_%d', i);
    model_name = sprintf('rem_%d', i);
    
    MODEL.Esc_rem.(fieldName).MODEL = M;
    MODEL.Esc_rem.(fieldName).name = model_name;
    
end

%Filtro de kalman para cada modelo

for i = 1:length(ss_rem_gdp)
    
    fieldName = sprintf('rem_%d', i);
    [MODEL.Esc_rem.(fieldName).MF, MODEL.Esc_rem.(fieldName).F] = filter(MODEL.Esc_rem.(fieldName).MODEL, MODEL.PreProc.obs,... 
                                MODEL.DATES.hist_start:MODEL.DATES.hist_end, ... 
                                'meanOnly=',true);
end

% Pronosticos para cada modelo

fcstrng = MODEL.DATES.pred_start:MODEL.DATES.pred_end;

for i = 1:length(ss_rem_gdp)
    
    fieldName = sprintf('rem_%d', i);
    MODEL.Esc_rem.(fieldName).F_pred = simulate(MODEL.Esc_rem.(fieldName).MF,... Modelo Filtrado
                            MODEL.Esc_rem.(fieldName).F,...Base de datos inicial filtrada
                            fcstrng,... Rango de Simulación
                            'anticipate', false,...
                            'DbOverlay=', true);

end
                    
for i = 1:length(ss_rem_gdp)
    
    fieldName = sprintf('rem_%d', i);
    plot(qq(2024,1):qq(2032,1),...
        MODEL.Esc_rem.(fieldName).F_pred.L_GDP_GAP);
    
    hold on
    
    zeroline();

end

MODEL.Esc_rem.rem_19.F_pred.L_GDP_GAP.plot
zeroline();


                    