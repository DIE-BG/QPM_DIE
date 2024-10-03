%% Escenario alterno 5
%{
    Escenarios de riesgo alrededor de variables externas o que no son
    instrumentos de Política Monetaria.

    Este Escenario: Para las variables mensuales se anclara el último valor
    observado mensual. 
%}

%% Anclaje de último dato mensual observado
MODEL.Esc.v5.name = MODEL.esc_names{6};
            
% Dato de anclaje
MODEL.DATES.E5_dates = MODEL.DATES.pred_start;
% return
%% %%%%%%%%%%%%%%%% Creación de escenario alternativo %%%%%%%%%%%%%%%%%%%%%%
% Shocks del escenario base
shocks = MODEL.F*get(MODEL.MF, 'elist');
% Concatenación de bd con condiciones iniciales (MODEL.F) y shocks Esc.
% Libre
MODEL.Esc.v5.dbi = dboverlay(MODEL.F,shocks);

% Anclaje de último valor observado mensual
% Inflación subyacente PCE de EEUU
MODEL.Esc.v5.dbi.D4L_CPI_RW(MODEL.DATES.E5_dates) = MODEL.PreProc.monthly.CPI_RW.log.diff(-12).data(end)*100;
%IPEI 
MODEL.Esc.v5.dbi.D4L_IPEI(MODEL.DATES.E5_dates) = MODEL.PreProc.monthly.IPEI.log.diff(-12).data(end)*100;
% Tasa de Fondos Federales
MODEL.Esc.v5.dbi.RS_RW(MODEL.DATES.E5_dates) = MODEL.PreProc.monthly.RS_RW(end);
% Inflación No Subyacente
MODEL.Esc.v5.dbi.D4L_CPI_NOSUBY(MODEL.DATES.E5_dates) = MODEL.PreProc.monthly.CPI.log.diff(-12).data(end)*100 - MODEL.PreProc.monthly.CPIXFE.log.diff(-12).data(end)*100;
% Inflación Subyacente
MODEL.Esc.v5.dbi.D4L_CPIXFE(MODEL.DATES.E5_dates) = MODEL.PreProc.monthly.CPIXFE.log.diff(-12).data(end)*100;
% Tipo de Cambio Nominal
MODEL.Esc.v5.dbi.D4L_S(MODEL.DATES.E5_dates) = MODEL.PreProc.monthly.S.log.diff(-12).data(end)*100;
% Base Monetaria
MODEL.Esc.v5.dbi.D4L_MB(MODEL.DATES.E5_dates) = MODEL.PreProc.monthly.MB.log.diff(-12).data(end)*100;
% Tasa de interés Lider
MODEL.Esc.v5.dbi.RS(MODEL.DATES.E5_dates) = MODEL.PreProc.monthly.RS(end);

%% Creación de plan de simulación
% Variables
vars = {'D4L_CPI_RW', 'D4L_IPEI', 'RS_RW','D4L_CPI_NOSUBY','D4L_CPIXFE','D4L_S','D4L_MB','RS'};
% choques
s_vars = {'SHK_DLA_CPI_RW', 'SHK_D4L_IPEI', 'SHK_RS_RW','SHK_D4L_CPI_NOSUBY','SHK_DLA_CPIXFE','SHK_L_S','SHK_D4L_MB','SHK_RS'};

% Plan de simulación
MODEL.Esc.v5.planSim = plan(MODEL.MF, MODEL.DATES.pred_start:MODEL.DATES.pred_end);
% Variable a endogenizar (shock propio?? No necesariamente)
MODEL.Esc.v5.planSim = endogenize(MODEL.Esc.v5.planSim, s_vars, MODEL.DATES.E5_dates); 
% Variable a exogenizar (Anclaje)
MODEL.Esc.v5.planSim = exogenize(MODEL.Esc.v5.planSim, vars, MODEL.DATES.E5_dates);

%% Simulación.
MODEL.Esc.v5.pred = simulate(MODEL.MF,...
                  MODEL.Esc.v5.dbi,...
                  MODEL.DATES.pred_start:MODEL.DATES.pred_end,...
                  'plan',MODEL.Esc.v5.planSim,...
                  'anticipate',false,...
                  'DbOverlay=', true);
% Descomposición             
MODEL.Esc.v5.shd = simulate(MODEL.MF,...
                  MODEL.Esc.v5.pred,...
                  MODEL.DATES.hist_start:MODEL.DATES.pred_end,...
                  'anticipate',false,...
                  'contributions',true);


%% Post-Procesamiento de variables seleccionadas.
% Desestacionalizar y obtener brechas y tendencias de estas variables
pp_list = {'L_MB', 'L_VEL', 'L_CPI_RW', 'L_CPI_RW_Q','L_Z', 'L_GDP', 'L_GDP_RW'};
% Recuperar niveles de estas variables
list_nivel = {'L_S','L_MB'};
                                        
MODEL = rec_GDP_RW(MODEL, 'Esc', 'v5');

MODEL = PostProcessing(MODEL,...
    'list',pp_list,...
    'list_niv', list_nivel,...
    'Esc',{'v5', MODEL.Esc.v5.pred});

%% Graficas

if graph_esc == true
    % Pre-processing
    % monthly
    PreProcPlots_m(MODEL,...
        'Esc_add', {'v5', MODEL_ANT},...
        'tab_range_mm', tab_range_mm);
    
    % quarterly
    PreProcPlots_q(MODEL,...
        'Esc_add', {'v5', MODEL_ANT},...
        'tab_range', tab_range_source_data)
    % Graficas de simulación
    simPlots(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v5','prediction_compared'),...}
        'Esc_add', {'v5', MODEL.Esc.v5.pred, MODEL.esc_col{5}},...
        'PlotList', get(MODEL.MF, 'xlist'),...
        'LegendsNames',{MODEL.esc_names{6}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Post Procesamiento
    % Logaritmo/tendencia vs corr. Anterior
    PostPrLogsComp(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE,'v5', 'PostProcessing'),...
        'Esc',{'v0', MODEL.PostProc},...
        'Esc_add', {'v5', MODEL.PostProc.v5, MODEL.esc_col{5}},...
        'PlotList', pp_list,...
        'LegendsNames',{MODEL.esc_names{6}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Variables en Niveles originales
    PostPrLevels(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v5','PostProcessing'),...
        'Esc_add', {'v5', MODEL.PostProc.v5, MODEL.esc_col{5}},...
        'PlotList', list_lev,...
        'Titles', tit_lev,...
        'LegendsNames',{MODEL.esc_names{6}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Brechas
    PostPrGaps(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v5','PostProcessing'),...
        'Esc_add', {'v5', MODEL.PostProc.v5, MODEL.esc_col{5}},...
        'PlotList', list_gaps,...
        'LegendsNames',{MODEL.esc_names{6}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Descomposición de choques para variables seleccionadas
    list = {'L_GDP_RW_GAP', 'DLA_CPI_RW', 'RS_RW', 'D4L_CPI_NOSUBY','L_GDP_GAP','DLA_CPIXFE', 'DLA_S', 'D4L_MB', 'RS',...
            'D4L_CPI', 'L_Z_GAP', 'D4L_VEL', 'RR', 'D4L_S', 'RMC', 'MCI'};
    
    % Long
    plot_shd_dsc(MODEL, MODEL.Esc.v5.pred, MODEL.Esc.v5.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v5','Shock_dec\long'),...
        'Rng', {},...
        'Variables',list)

    %short
    plot_shd_dsc(MODEL, MODEL.Esc.v5.pred, MODEL.Esc.v5.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v5','Shock_dec\short'),...
        'Rng', MODEL.DATES.hist_end-20:MODEL.DATES.hist_end+20,...
        'Variables',list)
    
    % Descomposición de chosques para variables seleccionadas
    % (Primeras diferencias)
    plot_diff_shd_dsc(MODEL,...
                  'variables', list,...
                  'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v5', 'Shock_dec', 'diff'),...
                  'Esc_add', {'v5', MODEL.Esc.v5.shd, MODEL.Esc.v5.pred});
    
    % Contribuciones
    contributions(MODEL,...
        'Esc_add', {'v5', MODEL.Esc.v5.pred}, ... 
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v5', 'contributions') ...           
    );
              
    contributions(MODEL,...
        'Esc_add', {'v5', MODEL.Esc.v5.pred}, ... 
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v5', 'diff_contributions'), ... 
        'Difference', true ...
    );             
    
    % Real exchange rate (subplot)
    tcr_subplot(MODEL,...
        'Esc_add', {'v5', MODEL.Esc.v5.pred, MODEL.esc_col{5}},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.esc_names{6}, MODEL.leg_act});
    
    % Money Velocity (subplot)
    vel_subplot(MODEL,...
        'tab_range', tab_range,...
        'Esc_add', {'v5', MODEL.PostProc, MODEL.esc_col{5}},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.esc_names{6}, MODEL.leg_act},...
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v5', 'otras'));
end
disp('Escenario 5: ok');