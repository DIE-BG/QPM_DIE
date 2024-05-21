%% Escenario Contrafactual de PM
%{
    Escenarios contrafactuales utilizando instrumentos de PM:
                - Tasa de interes líder
                - Base Monetaria

    No necesariamente se "endogeniza" el choque correspondiente a la
    variable anclada, ni en los mismos períodos.
%}

%% 1. Tasa de interés invariable en el primer trimestre de pronóstico

MODEL.Esc.v2.dbi = MODEL.F;

% Mantener la tasa invariable
MODEL.Esc.v2.dbi.RS(MODEL.DATES.pred_start) = MODEL.Esc.v2.dbi.RS(MODEL.DATES.hist_end);

% Plan de simulación
MODEL.Esc.v2.planSim = plan(MODEL.MF,MODEL.DATES.pred_start:MODEL.DATES.pred_end);
% Variable a endogenizar (shock propio?? No necesariamente)
MODEL.Esc.v2.planSim = endogenize(MODEL.Esc.v2.planSim,{'SHK_RS'},MODEL.DATES.pred_start); 
% Variable a exogenizar (Anclaje)
MODEL.Esc.v2.planSim = exogenize(MODEL.Esc.v2.planSim,{'RS'},MODEL.DATES.pred_start);

%% Simulación.
MODEL.Esc.v2.pred = simulate(MODEL.MF,...
              MODEL.Esc.v2.dbi,...
              MODEL.DATES.pred_start:MODEL.DATES.pred_end,...
              'plan',MODEL.Esc.v2.planSim,...
              'anticipate',false,...
              'DbOverlay=', true);     
          
%% Descomposición             
MODEL.Esc.v2.shd = simulate(MODEL.MF,...
                  MODEL.Esc.v2.pred,...
                  MODEL.DATES.hist_start:MODEL.DATES.pred_end,...
                  'anticipate',false,...
                  'contributions',true);
              
%% Post-Procesamiento de variables seleccionadas.
pp_list = {'L_MB', 'L_VEL', 'L_CPI_RW', 'L_CPI_RW_Q','L_Z', 'L_GDP', 'L_GDP_RW'};
list_nivel = {'L_S','L_MB'};

MODEL = rec_GDP_RW(MODEL, 'Esc', 'v2');
                                        
MODEL = PostProcessing(MODEL,...
    'list',pp_list,...
    'list_niv', list_nivel,...
    'Esc',{'v2', MODEL.Esc.v2.pred});            

%% Graficas

    % Pre-processing
    % monthly
    PreProcPlots_m(MODEL,...
        'Esc_add', {'v2', MODEL_ANT},...
        'tab_range_mm', tab_range_mm);
    
    % quarterly
    PreProcPlots_q(MODEL,...
        'Esc_add', {'v2', MODEL_ANT},...
        'tab_range', tab_range_source_data)
    % Graficas de simulación
    simPlots(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v2','prediction_compared'),...}
        'Esc_add', {'v2', MODEL.Esc.v2.pred, MODEL.esc_col{1}},...
        'PlotList', get(MODEL.MF, 'xlist'),...
        'LegendsNames',{MODEL.esc_names{3}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Post Procesamiento
    % Logaritmo/tendencia vs corr. Anterior
    PostPrLogsComp(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE,'v2', 'PostProcessing'),...
        'Esc',{'v0', MODEL.PostProc},...
        'Esc_add', {'v2', MODEL.PostProc.v2, MODEL.esc_col{1}},...
        'PlotList', pp_list,...
        'LegendsNames',{MODEL.esc_names{3}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Variables en Niveles originales
    PostPrLevels(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v2','PostProcessing'),...
        'Esc_add', {'v2', MODEL.PostProc.v2, MODEL.esc_col{1}},...
        'PlotList', list_lev,...
        'Titles', tit_lev,...
        'LegendsNames',{MODEL.esc_names{3}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Brechas
    PostPrGaps(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v2','PostProcessing'),...
        'Esc_add', {'v2', MODEL.PostProc.v2, MODEL.esc_col{1}},...
        'PlotList', list_gaps,...
        'LegendsNames',{MODEL.esc_names{3}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Descomposición de choques para variables seleccionadas
    list = {'L_GDP_RW_GAP', 'DLA_CPI_RW', 'RS_RW', 'D4L_CPI_NOSUBY','L_GDP_GAP','DLA_CPIXFE', 'DLA_S', 'D4L_MB', 'RS',...
            'D4L_CPI', 'L_Z_GAP', 'D4L_VEL', 'RR', 'D4L_S', 'RMC', 'MCI'};
    
    % Long
    plot_shd_dsc(MODEL, MODEL.Esc.v2.pred, MODEL.Esc.v2.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v2','Shock_dec\long'),...
        'Rng', {},...
        'Variables',list)

    %short
    plot_shd_dsc(MODEL, MODEL.Esc.v2.pred, MODEL.Esc.v2.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v2','Shock_dec\short'),...
        'Rng', MODEL.DATES.hist_end-20:MODEL.DATES.hist_end+20,...
        'Variables',list)
    
    % Descomposición de chosques para variables seleccionadas
    % (Primeras diferencias)
    plot_diff_shd_dsc(MODEL,...
                  'variables', list,...
                  'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v2', 'Shock_dec', 'diff'),...
                  'Esc_add', {'v2', MODEL.Esc.v2.shd, MODEL.Esc.v2.pred});
    
    % Contribuciones
    contributions(MODEL,...
                  'Esc_add', {'v2', MODEL.Esc.v2.pred});
   diff_contributions(MODEL,...
                  'Esc_add', {'v2', MODEL.Esc.v2.pred}); 
    
    % Real exchange rate (subplot)
    tcr_subplot(MODEL,...
        'Esc_add', {'v2', MODEL.Esc.v2.pred, MODEL.esc_col{1}},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.esc_names{3}, MODEL.leg_act});
    
    % Money Velocity (subplot)
    vel_subplot(MODEL,...
        'tab_range', tab_range,...
        'Esc_add', {'v2', MODEL.PostProc, MODEL.esc_col{1}},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.esc_names{3}, MODEL.leg_act},...
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v2', 'otras'));
    
disp('Escenario 2: ok');