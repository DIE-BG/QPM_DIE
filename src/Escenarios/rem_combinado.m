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

MODEL.Esc.v6.name = 'Escenario remesas combinado';
MODEL.Esc.v6.pred = MODEL.F_REM_pred;
MODEL.Esc.v6.dbi = MODEL.F_REM;

% Descopmposición de shocks
MODEL.Esc.v6.shd = simulate(MODEL.MF_REM,...
                  MODEL.Esc.v6.pred,...
                  MODEL.DATES.hist_start:MODEL.DATES.pred_end,...
                  'anticipate',false,...
                  'contributions',true);
                    


%% Post-Procesamiento de variables seleccionadas.
% Desestacionalizar y obtener brechas y tendencias de estas variables
pp_list = {'L_MB', 'L_VEL', 'L_CPI_RW', 'L_IPEI_Q','L_Z', 'L_GDP', 'L_GDP_RW'};
% Recuperar niveles de estas variables
list_nivel = {'L_S','L_MB'};
                                        
MODEL = rec_GDP_RW(MODEL, 'Esc', 'v6');

MODEL = PostProcessing(MODEL,...
    'list',pp_list,...
    'list_niv', list_nivel,...
    'Esc',{'v6', MODEL.Esc.v6.pred});

%% Graficas del escenario combinado
col = [0.89 0.212 0.843];
graph_esc = true;
if graph_esc == true
    % Pre-processing
    % monthly
    PreProcPlots_m(MODEL,...
        'Esc_add', {'v6', MODEL_ANT},...
        'tab_range_mm', tab_range_mm);
    
    % quarterly
    PreProcPlots_q(MODEL,...
        'Esc_add', {'v6', MODEL_ANT},...
        'tab_range', tab_range_source_data)
    % Graficas de simulación
    simPlots(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v6','prediction_compared'),...}
        'Esc_add', {'v6', MODEL.Esc.v6.pred, col},...
        'PlotList', get(MODEL.MF, 'xlist'),...
        'LegendsNames',{MODEL.Esc.v6.name, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Post Procesamiento
    % Logaritmo/tendencia vs corr. Anterior
    PostPrLogsComp(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v6','PostProcessing'),...
        'Esc',{'v0', MODEL.PostProc},...
        'Esc_add', {'v6', MODEL.PostProc.v6, col},...
        'PlotList', pp_list,...
        'LegendsNames',{MODEL.Esc.v6.name, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Variables en Niveles originales
    PostPrLevels(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v6','PostProcessing'),...
        'Esc_add', {'v6', MODEL.PostProc.v6, col},...
        'PlotList', list_lev,...
        'Titles', tit_lev,...
        'LegendsNames',{MODEL.Esc.v6.name, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Brechas
    PostPrGaps(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v6','PostProcessing'),...
        'Esc_add', {'v6', MODEL.PostProc.v6, col},...
        'PlotList', list_gaps,...
        'LegendsNames',{MODEL.Esc.v6.name, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Descomposición de choques para variables seleccionadas
    list = {'L_GDP_RW_GAP', 'DLA_CPI_RW', 'RS_RW', 'D4L_CPI_NOSUBY','L_GDP_GAP','DLA_CPIXFE', 'DLA_S', 'D4L_MB', 'RS',...
            'D4L_CPI', 'L_Z_GAP', 'D4L_VEL', 'RR', 'D4L_S', 'RMC', 'MCI'};
    
    % Long
    plot_shd_dsc(MODEL, MODEL.Esc.v6.pred, MODEL.Esc.v6.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v6','Shock_dec\long'),...
        'Rng', {},...
        'Variables',list)

    %short
    plot_shd_dsc(MODEL, MODEL.Esc.v6.pred, MODEL.Esc.v6.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v6','Shock_dec\short'),...
        'Rng', MODEL.DATES.hist_end-20:MODEL.DATES.hist_end+20,...
        'Variables',list)
    
    % Descomposición de chosques para variables seleccionadas
    % (Primeras diferencias)
    plot_diff_shd_dsc(MODEL,...
                  'variables', list,...
                  'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v6', 'Shock_dec', 'diff'),...
                  'Esc_add', {'v6', MODEL.Esc.v6.shd, MODEL.Esc.v6.pred});
    
    % Contribuciones
    contributions(MODEL,...
        'Esc_add', {'v6', MODEL.Esc.v6.pred}, ... 
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v6', 'contributions') ...           
    );
              
    contributions(MODEL,...
        'Esc_add', {'v6', MODEL.Esc.v6.pred}, ... 
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v6', 'diff_contributions'), ... 
        'Difference', true ...
    );             
    
    % Real exchange rate (subplot)
    tcr_subplot(MODEL,...
        'Esc_add', {'v6', MODEL.Esc.v6.pred, col},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.Esc.v6.name, MODEL.leg_act});
    
    % Money Velocity (subplot)
    vel_subplot(MODEL,...
        'tab_range', tab_range,...
        'Esc_add', {'v6', MODEL.PostProc, col},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.Esc.v6.name, MODEL.leg_act},...
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v6', 'otras'));
end
disp('Escenario 6: ok');



