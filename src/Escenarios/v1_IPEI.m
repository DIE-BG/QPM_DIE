%% Escenario alterno 1
%{
    Escenarios de riesgo alrededor de variables externas o que no son
    instrumentos de Política Monetaria.

    Este Escenario: Se anclan pronósticos de IPEI. 
    Se usan los pronósticos del logaritmo desestacionalizado provenientes 
    del post-procesamiento del SVAR50-4b
%}

%% Anclaje de IPEI proveniente de SVAR en horizonte de pronóstico
% (8 Trimestres)

alt1 = load(fullfile('data', 'corrimientos',MODEL.CORR_DATE,...
            'v1', sprintf('MODEL-%s.mat',MODEL.CORR_DATE)));

alt1 = alt1.MODEL.F_pred;

% Trimestres de anclaje
MODEL.DATES.E1_dates = MODEL.DATES.pred_start:MODEL.DATES.pred_start+7;



%% %%%%%%%%%%%%%%%% Creación de escenario alternativo %%%%%%%%%%%%%%%%%%%%%%
% Shocks del escenario base
shocks = MODEL.F*get(MODEL.MF, 'elist');
% Concatenación de bd con condiciones iniciales (MODEL.F) y shocks Esc.
% Libre
MODEL.Esc.v1.dbi = dboverlay(MODEL.F,shocks);

% Imposición de anclajes provenientes del QPM en base de datos
MODEL.Esc.v1.dbi.D4L_IPEI(MODEL.DATES.E1_dates) = alt1.d4_ln_ipei(MODEL.DATES.E1_dates);

% Plan de simulación
MODEL.Esc.v1.planSim = plan(MODEL.MF, MODEL.DATES.pred_start:MODEL.DATES.pred_end);
% Variable a endogenizar (shock propio?? No necesariamente)
MODEL.Esc.v1.planSim = endogenize(MODEL.Esc.v1.planSim,{'SHK_D4L_IPEI'}, MODEL.DATES.E1_dates); 
% Variable a exogenizar (Anclaje)
MODEL.Esc.v1.planSim = exogenize(MODEL.Esc.v1.planSim,{'D4L_IPEI'}, MODEL.DATES.E1_dates);

% Simulación.
MODEL.Esc.v1.pred = simulate(MODEL.MF,...
                  MODEL.Esc.v1.dbi,...
                  MODEL.DATES.pred_start:MODEL.DATES.pred_end,...
                  'plan',MODEL.Esc.v1.planSim,...
                  'anticipate',false,...
                  'DbOverlay=', true);

% Descomposición             
MODEL.Esc.v1.shd = simulate(MODEL.MF,...
                  MODEL.Esc.v1.pred,...
                  MODEL.DATES.hist_start:MODEL.DATES.pred_end,...
                  'anticipate',false,...
                  'contributions',true);



%% Post-Procesamiento de variables seleccionadas.
% Desestacionalizar y obtener brechas y tendencias de estas variables
pp_list = {'L_MB', 'L_VEL', 'L_CPI_RW', 'L_IPEI_Q','L_Z', 'L_GDP', 'L_GDP_RW'};
% Recuperar niveles de estas variables
list_nivel = {'L_S','L_MB'};
                                        
MODEL = rec_GDP_RW(MODEL, 'Esc', 'v1');

MODEL = PostProcessing(MODEL,...
    'list',pp_list,...
    'list_niv', list_nivel,...
    'Esc',{'v1', MODEL.Esc.v1.pred});

%% Graficas

if graph_esc == true
    % Pre-processing
    % monthly
    PreProcPlots_m(MODEL,...
        'Esc_add', {'v1', MODEL_ANT},...
        'tab_range_mm', tab_range_mm);
    
    % quarterly
    PreProcPlots_q(MODEL,...
        'Esc_add', {'v1', MODEL_ANT},...
        'tab_range', tab_range_source_data)
    % Graficas de simulación
    simPlots(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v1','prediction_compared'),...}
        'Esc_add', {'v1', MODEL.Esc.v1.pred, MODEL.esc_col{2}},...
        'PlotList', get(MODEL.MF, 'xlist'),...
        'LegendsNames',{MODEL.esc_names{2}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Post Procesamiento
    % Logaritmo/tendencia vs corr. Anterior
    PostPrLogsComp(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE,'v1', 'PostProcessing'),...
        'Esc',{'v0', MODEL.PostProc},...
        'Esc_add', {'v1', MODEL.PostProc.v1, MODEL.esc_col{2}},...
        'PlotList', pp_list,...
        'LegendsNames',{MODEL.esc_names{2}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Variables en Niveles originales
    PostPrLevels(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v1','PostProcessing'),...
        'Esc_add', {'v1', MODEL.PostProc.v1, MODEL.esc_col{2}},...
        'PlotList', list_lev,...
        'Titles', tit_lev,...
        'LegendsNames',{MODEL.esc_names{2}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Brechas
    PostPrGaps(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v1','PostProcessing'),...
        'Esc_add', {'v1', MODEL.PostProc.v1, MODEL.esc_col{2}},...
        'PlotList', list_gaps,...
        'LegendsNames',{MODEL.esc_names{2}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Descomposición de choques para variables seleccionadas
    list = {'L_GDP_RW_GAP', 'DLA_CPI_RW', 'RS_RW', 'D4L_CPI_NOSUBY','L_GDP_GAP','DLA_CPIXFE', 'DLA_S', 'D4L_MB', 'RS',...
            'D4L_CPI', 'L_Z_GAP', 'D4L_VEL', 'RR', 'D4L_S', 'RMC', 'MCI'};
    
    % Long
    plot_shd_dsc(MODEL, MODEL.Esc.v1.pred, MODEL.Esc.v1.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v1','Shock_dec\long'),...
        'Rng', {},...
        'Variables',list)

    %short
    plot_shd_dsc(MODEL, MODEL.Esc.v1.pred, MODEL.Esc.v1.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v1','Shock_dec\short'),...
        'Rng', MODEL.DATES.hist_end-20:MODEL.DATES.hist_end+20,...
        'Variables',list)
    
    % Descomposición de chosques para variables seleccionadas
    % (Primeras diferencias)
    plot_diff_shd_dsc(MODEL,...
                  'variables', list,...
                  'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v1', 'Shock_dec', 'diff'),...
                  'Esc_add', {'v1', MODEL.Esc.v1.shd, MODEL.Esc.v1.pred});
    
    % Contribuciones
    contributions(MODEL,...
        'Esc_add', {'v1', MODEL.Esc.v1.pred}, ... 
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v1', 'contributions') ...           
    );
              
    contributions(MODEL,...
        'Esc_add', {'v1', MODEL.Esc.v1.pred}, ... 
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v1', 'diff_contributions'), ... 
        'Difference', true ...
    );         
    
    % Real exchange rate (subplot)
    tcr_subplot(MODEL,...
        'Esc_add', {'v1', MODEL.Esc.v1.pred, MODEL.esc_col{2}},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.esc_names{2}, MODEL.leg_act});
    
    % Money Velocity (subplot)
    vel_subplot(MODEL,...
        'tab_range', tab_range,...
        'Esc_add', {'v1', MODEL.PostProc, MODEL.esc_col{2}},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.esc_names{2}, MODEL.leg_act},...
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v1', 'otras'));
end
disp('Escenario 1: ok');