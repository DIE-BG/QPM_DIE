
%% Anclaje de IPEI proveniente de SVAR en horizonte de pronóstico
% (8 Trimestres)
MODEL.Esc.v3.name = MODEL.esc_names{4};
alt3 = load(fullfile('data', 'corrimientos',MODEL.CORR_DATE,...
            'v1', sprintf('MODEL-%s.mat',MODEL.CORR_DATE)));

alt3 = alt3.MODEL.PostProc.v0.l_sa;
        
% Trimestres de anclaje
MODEL.DATES.E3_dates = MODEL.DATES.E1_dates;

%% %%%%%%%%%%%%%%%% Creación de escenario alternativo %%%%%%%%%%%%%%%%%%%%%%
% Shocks del escenario base
shocks = MODEL.F*get(MODEL.MF, 'elist');
MODEL.Esc.v3.dbi = dboverlay(MODEL.F,shocks);

% Imposición de anclajes provenientes del QPM en base de datos
MODEL.Esc.v3.dbi.L_CPI_RW(MODEL.DATES.E3_dates) = alt3.ln_ipei_sa(MODEL.DATES.E3_dates);
% Mantenimiento de Tasa un trimestre
MODEL.Esc.v3.dbi.RS(MODEL.DATES.pred_start) = MODEL.Esc.v3.dbi.RS(MODEL.DATES.hist_end);

% Plan de simulación
MODEL.Esc.v3.planSim = plan(MODEL.MF, MODEL.DATES.pred_start:MODEL.DATES.pred_end);
% Variable a endogenizar (shock propio?? No necesariamente)
MODEL.Esc.v3.planSim = endogenize(MODEL.Esc.v3.planSim,{'SHK_DLA_CPI_RW'},MODEL.DATES.E1_dates); 
MODEL.Esc.v3.planSim = endogenize(MODEL.Esc.v3.planSim,{'SHK_RS'},MODEL.DATES.pred_start); 
% Variable a exogenizar (Anclaje)
MODEL.Esc.v3.planSim = exogenize(MODEL.Esc.v3.planSim,{'L_CPI_RW'},MODEL.DATES.E1_dates);
MODEL.Esc.v3.planSim = exogenize(MODEL.Esc.v3.planSim,{'RS'},MODEL.DATES.pred_start);

%% Simulación.
MODEL.Esc.v3.pred = simulate(MODEL.MF,...
                  MODEL.Esc.v3.dbi,...
                  MODEL.DATES.pred_start:MODEL.DATES.pred_end,...
                  'plan',MODEL.Esc.v3.planSim,...
                  'anticipate',false,...
                  'DbOverlay=', true);

% Descomposición             
MODEL.Esc.v3.shd = simulate(MODEL.MF,...
                  MODEL.Esc.v3.pred,...
                  MODEL.DATES.hist_start:MODEL.DATES.pred_end,...
                  'anticipate',false,...
                  'contributions',true);



%% Post-Procesamiento de variables seleccionadas.
pp_list = {'L_MB', 'L_VEL', 'L_CPI_RW', 'L_CPI_RW_Q','L_Z', 'L_GDP', 'L_GDP_RW'};
list_nivel = {'L_S','L_MB'};
                                        
MODEL = rec_GDP_RW(MODEL, 'Esc', 'v3');

MODEL = PostProcessing(MODEL,...
    'list',pp_list,...
    'list_niv', list_nivel,...
    'Esc',{'v3', MODEL.Esc.v3.pred});


%% Graficas
if graph_esc == true
    % Pre-processing
    % monthly
    PreProcPlots_m(MODEL,...
        'Esc_add', {'v3', MODEL_ANT},...
        'tab_range_mm', tab_range_mm);
    
    % quarterly
    PreProcPlots_q(MODEL,...
        'Esc_add', {'v3', MODEL_ANT},...
        'tab_range', tab_range_source_data)
    % Graficas de simulación
    simPlots(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v3','prediction_compared'),...}
        'Esc_add', {'v3', MODEL.Esc.v3.pred, MODEL.esc_col{3}},...
        'PlotList', get(MODEL.MF, 'xlist'),...
        'LegendsNames',{MODEL.esc_names{4}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Post Procesamiento
    % Logaritmo/tendencia vs corr. Anterior
    PostPrLogsComp(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE,'v3', 'PostProcessing'),...
        'Esc',{'v0', MODEL.PostProc},...
        'Esc_add', {'v3', MODEL.PostProc.v3, MODEL.esc_col{3}},...
        'PlotList', pp_list,...
        'LegendsNames',{MODEL.esc_names{4}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Variables en Niveles originales
    PostPrLevels(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v3','PostProcessing'),...
        'Esc_add', {'v3', MODEL.PostProc.v3, MODEL.esc_col{3}},...
        'PlotList', list_lev,...
        'Titles', tit_lev,...
        'LegendsNames',{MODEL.esc_names{4}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Brechas
    PostPrGaps(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v3','PostProcessing'),...
        'Esc_add', {'v3', MODEL.PostProc.v3, MODEL.esc_col{3}},...
        'PlotList', list_gaps,...
        'LegendsNames',{MODEL.esc_names{4}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Descomposición de choques para variables seleccionadas
    list = {'L_GDP_RW_GAP', 'DLA_CPI_RW', 'RS_RW', 'D4L_CPI_NOSUBY','L_GDP_GAP','DLA_CPIXFE', 'DLA_S', 'D4L_MB', 'RS',...
            'D4L_CPI', 'L_Z_GAP', 'D4L_VEL', 'RR', 'D4L_S', 'RMC', 'MCI'};
    
    % Long
    plot_shd_dsc(MODEL, MODEL.Esc.v3.pred, MODEL.Esc.v3.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v3','Shock_dec\long'),...
        'Rng', {},...
        'Variables',list)

    %short
    plot_shd_dsc(MODEL, MODEL.Esc.v3.pred, MODEL.Esc.v3.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v3','Shock_dec\short'),...
        'Rng', MODEL.DATES.hist_end-20:MODEL.DATES.hist_end+20,...
        'Variables',list)
    
    % Descomposición de chosques para variables seleccionadas
    % (Primeras diferencias)
    plot_diff_shd_dsc(MODEL,...
                  'variables', list,...
                  'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v3', 'Shock_dec', 'diff'),...
                  'Esc_add', {'v3', MODEL.Esc.v3.shd, MODEL.Esc.v3.pred});
    
    % Contribuciones
    contributions(MODEL,...
                  'Esc_add', {'v3', MODEL.Esc.v3.pred});
    diff_contributions(MODEL,...
                  'Esc_add', {'v3', MODEL.Esc.v3.pred});     
    % Real exchange rate (subplot)
    tcr_subplot(MODEL,...
        'Esc_add', {'v3', MODEL.Esc.v3.pred, MODEL.esc_col{3}},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.esc_names{4}, MODEL.leg_act});
    
    % Money Velocity (subplot)
    vel_subplot(MODEL,...
        'tab_range', tab_range,...
        'Esc_add', {'v3', MODEL.PostProc, MODEL.esc_col{3}},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.esc_names{4}, MODEL.leg_act},...
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v3', 'otras'));
    
    
        MODEL = fanchart(MODEL,'PlotList',{'D4L_CPI','L_GDP_GAP','RS'},...
                    'EndDatePlot',qq(2025,4),...
                    'Esc', 'v3',...
                    'apertura', {[0.25 0.50 0.75 0.75 0.75 0.75 0.75]',...
                                {},{}},... Vector Columna
                    'Grilla',[0.05:0.05:0.95],...
                    'sesgo', {[1 1 1 1 1 1 1 1 1,... Percentiles por debajo de la mediana
                               1,... percentil 50
                               1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5],... percentiles por arriba de la mediana
                               {},{}}...Vector fila
                    );      
end
disp('Escenario 23: ok');