%{
Compartativo moficiando el QPM_MA añadiendo (1-RHO_REM_GDP)*REM_GDP_BAR a
la ecuación de REM_GDP
%}

% parámetros del QPM
setparam;

% Modelo
M = model('QPM_MA_modif.model', 'assign', s);
M = sstate(M,'growth=',true,'MaxFunEvals',1000,'display=','off');
M = solve(M,'error=',true);

MODEL.Esc.v8.MODEL = M;

% Filtro de kalman
[MODEL.Esc.v8.MF, MODEL.Esc.v8.F] = filter(MODEL.Esc.v8.MODEL, MODEL.PreProc.obs,... 
                                MODEL.DATES.hist_start:MODEL.DATES.hist_end, ... 
                                'meanOnly=',true);
                            
% Pronosticos
fcstrng = MODEL.DATES.pred_start:MODEL.DATES.pred_end;
MODEL.Esc.v8.pred = simulate(MODEL.Esc.v8.MF,... Modelo Filtrado
                        MODEL.Esc.v8.F,...Base de datos inicial filtrada
                        fcstrng,... Rango de Simulación
                        'anticipate', false,...
                        'DbOverlay=', true);

% Descopmposición de shocks
MODEL.Esc.v8.shd = simulate(MODEL.Esc.v8.MF,...
                  MODEL.Esc.v8.pred,...
                  MODEL.DATES.hist_start:MODEL.DATES.pred_end,...
                  'anticipate',false,...
                  'contributions',true);

%% Post-Procesamiento de variables seleccionadas.
% Desestacionalizar y obtener brechas y tendencias de estas variables
pp_list = {'L_MB', 'L_VEL', 'L_CPI_RW', 'L_IPEI_Q','L_Z', 'L_GDP', 'L_GDP_RW'};
% Recuperar niveles de estas variables
list_nivel = {'L_S','L_MB'};
                                        
MODEL = rec_GDP_RW(MODEL, 'Esc', 'v8');

MODEL = PostProcessing(MODEL,...
    'list',pp_list,...
    'list_niv', list_nivel,...
    'Esc',{'v8', MODEL.Esc.v8.pred});

%% Cambio de etiqueta para el modelo con media movil
leg_act = 'Octubre 2024 rem_gdp_bar MA'; 
MODEL.Esc.v8.name = 'Octubre 2024 rem_gdp_bar MA (modificado)';

%% Graficas del escenario combinado
col = [0.89 0.212 0.843];
graph_esc = true;
if graph_esc == true
    % Pre-processing
    % monthly
    PreProcPlots_m(MODEL,...
        'Esc_add', {'v8', MODEL_ANT},...
        'tab_range_mm', tab_range_mm);
    
    % quarterly
    PreProcPlots_q(MODEL,...
        'Esc_add', {'v8', MODEL_ANT},...
        'tab_range', tab_range_source_data)
    % Graficas de simulación
    simPlots(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v8','prediction_compared'),...}
        'Esc_add', {'v8', MODEL.Esc.v8.pred, col},...
        'PlotList', get(MODEL.MF, 'xlist'),...
        'LegendsNames',{MODEL.Esc.v8.name, leg_act},...
        'TabRange', tab_range...
        );
    
    % Post Procesamiento
    % Logaritmo/tendencia vs corr. Anterior
    PostPrLogsComp(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v8','PostProcessing'),...
        'Esc',{'v0', MODEL.PostProc},...
        'Esc_add', {'v8', MODEL.PostProc.v8, col},...
        'PlotList', pp_list,...
        'LegendsNames',{MODEL.Esc.v8.name, leg_act},...
        'TabRange', tab_range...
        );
    % Variables en Niveles originales
    PostPrLevels(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v8','PostProcessing'),...
        'Esc_add', {'v8', MODEL.PostProc.v8, col},...
        'PlotList', list_lev,...
        'Titles', tit_lev,...
        'LegendsNames',{MODEL.Esc.v8.name, leg_act},...
        'TabRange', tab_range...
        );
    % Brechas
    PostPrGaps(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v8','PostProcessing'),...
        'Esc_add', {'v8', MODEL.PostProc.v8, col},...
        'PlotList', list_gaps,...
        'LegendsNames',{MODEL.Esc.v8.name, leg_act},...
        'TabRange', tab_range...
        );
    
    % Descomposición de choques para variables seleccionadas
    list = {'L_GDP_RW_GAP', 'DLA_CPI_RW', 'RS_RW', 'D4L_CPI_NOSUBY','L_GDP_GAP','DLA_CPIXFE', 'DLA_S', 'D4L_MB', 'RS',...
            'D4L_CPI', 'L_Z_GAP', 'D4L_VEL', 'RR', 'D4L_S', 'RMC', 'MCI'};
    
    % Long
    plot_shd_dsc(MODEL, MODEL.Esc.v8.pred, MODEL.Esc.v8.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v8','Shock_dec\long'),...
        'MF', MODEL.Esc.v8.MF, 'Rng', {},...
        'Variables',list)

    %short
    plot_shd_dsc(MODEL, MODEL.Esc.v8.pred, MODEL.Esc.v8.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v8','Shock_dec\short'),...
        'MF', MODEL.Esc.v8.MF, 'Rng', MODEL.DATES.hist_end-20:MODEL.DATES.hist_end+20,...
        'Variables',list)
    
    % Descomposición de chosques para variables seleccionadas
    % (Primeras diferencias)
    plot_diff_shd_dsc(MODEL,...
                  'variables', list,...
                  'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v8', 'Shock_dec', 'diff'),...
                  'MF', MODEL.Esc.v8.MF, 'Esc_add', {'v8', MODEL.Esc.v8.shd, MODEL.Esc.v8.pred});
    
    % Contribuciones
    contributions(MODEL,...
        'Esc_add', {'v8', MODEL.Esc.v8.pred}, ... 
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v8', 'contributions') ...           
    );
              
    contributions(MODEL,...
        'Esc_add', {'v8', MODEL.Esc.v8.pred}, ... 
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v8', 'diff_contributions'), ... 
        'Difference', true ...
    );             
    
    % Real exchange rate (subplot)
    tcr_subplot(MODEL,...
        'Esc_add', {'v8', MODEL.Esc.v8.pred, col},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.Esc.v8.name, leg_act});
    
    % Money Velocity (subplot)
    vel_subplot(MODEL,...
        'tab_range', tab_range,...
        'Esc_add', {'v8', MODEL.PostProc, col},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.Esc.v8.name, leg_act},...
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v8', 'otras'));
end
disp('Escenario 8: ok');              