%% Escenario alterno 7
%{
    Escenarios de riesgo alrededor de variables externas o que no son
    instrumentos de Política Monetaria.

    Este Escenario: Se anclan pronósticos de la tasa de interes de fondos federales. 
    Se usan los pronósticos  provenientes del FOMC (se anclan los
    pronósticos medios de la Federal funds rate)
%}

%% Anclaje de tasa de interes de fondos federales y PCE CORE
% (Anclaje hasta 2025Q4)

% Base de datos
alt7 = databank.fromCSV(...
        fullfile('data', 'corrimientos',MODEL.CORR_DATE,...
                 'v6', 'fulldata_EXT.csv')); 
             
% Trimestres de anclaje
MODEL.DATES.E7_dates = MODEL.DATES.pred_start:MODEL.DATES.pred_start+5;
MODEL.DATES.E8_dates = MODEL.DATES.pred_start;

%% %%%%%%%%%%%%%%%% Creación de escenario alternativo %%%%%%%%%%%%%%%%%%%%%%
% Shocks del escenario base
shocks = MODEL.F*get(MODEL.MF, 'elist');
% Concatenación de bd con condiciones iniciales (MODEL.F) y shocks Esc.
% Libre
MODEL.Esc.v7.dbi = dboverlay(MODEL.F,shocks);

% Anclaje de último valor observado mensual
% Tipo de Cambio Nominal
MODEL.Esc.v7.dbi.D4L_S(MODEL.DATES.E8_dates) = MODEL.PreProc.monthly.S.log.diff(-12).data(end)*100;

% Imposición de anclajes provenientes del FOMC últimos pronósticos
% publicados
MODEL.Esc.v7.dbi.RS_RW(MODEL.DATES.E7_dates) = alt7.istar(MODEL.DATES.E7_dates);
MODEL.Esc.v7.dbi.DLA_CPI_RW(MODEL.DATES.E7_dates) = alt7.pceus(MODEL.DATES.E7_dates);
% MODEL.Esc.v7.dbi.DLA_GDP_BAR(MODEL.DATES.E7_dates) = alt7.gdppot(MODEL.DATES.E7_dates);

%% Creación de plan de simulación
% Variables
vars = {'RS_RW','DLA_CPI_RW'};
vars_2 = {'D4L_S'};
% choques
s_vars = {'SHK_RS_RW','SHK_DLA_CPI_RW'};
s_vars_2 =  {'SHK_L_S'};
% Plan de simulación
MODEL.Esc.v7.planSim = plan(MODEL.MF, MODEL.DATES.pred_start:MODEL.DATES.pred_end);
% Variable a endogenizar (shock propio?? No necesariamente)
MODEL.Esc.v7.planSim = endogenize(MODEL.Esc.v7.planSim,s_vars,MODEL.DATES.E7_dates); 
MODEL.Esc.v7.planSim = endogenize(MODEL.Esc.v7.planSim,s_vars_2,MODEL.DATES.E8_dates); 
% Variable a exogenizar (Anclaje)
MODEL.Esc.v7.planSim = exogenize(MODEL.Esc.v7.planSim,vars,MODEL.DATES.E7_dates);
MODEL.Esc.v7.planSim = exogenize(MODEL.Esc.v7.planSim,vars_2,MODEL.DATES.E8_dates);

%% Simulación.
MODEL.Esc.v7.pred = simulate(MODEL.MF,...
                  MODEL.Esc.v7.dbi,...
                  MODEL.DATES.pred_start:MODEL.DATES.pred_end,...
                  'plan',MODEL.Esc.v7.planSim,...
                  'anticipate',false,...
                  'DbOverlay=', true);

% Descomposición             
MODEL.Esc.v7.shd = simulate(MODEL.MF,...
                  MODEL.Esc.v7.pred,...
                  MODEL.DATES.hist_start:MODEL.DATES.pred_end,...
                  'anticipate',false,...
                  'contributions',true);



%% Post-Procesamiento de variables seleccionadas.
pp_list = {'L_MB', 'L_VEL', 'L_CPI_RW', 'L_IPEI_Q','L_Z', 'L_GDP', 'L_GDP_RW'};
list_nivel = {'L_S','L_MB'};
                                        
MODEL = rec_GDP_RW(MODEL, 'Esc', 'v7');

MODEL = PostProcessing(MODEL,...
    'list',pp_list,...
    'list_niv', list_nivel,...
    'Esc',{'v7', MODEL.Esc.v7.pred});

%% Graficas

if graph_esc == true
    % Pre-processing
    % monthly
    PreProcPlots_m(MODEL,...
        'Esc_add', {'v7', MODEL_ANT},...
        'tab_range_mm', tab_range_mm);
    
    % quarterly
    PreProcPlots_q(MODEL,...
        'Esc_add', {'v7', MODEL_ANT},...
        'tab_range', tab_range_source_data)
    % Graficas de simulación
    simPlots(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v7','prediction_compared'),...}
        'Esc_add', {'v7', MODEL.Esc.v7.pred, MODEL.esc_col{7}},...
        'PlotList', get(MODEL.MF, 'xlist'),...
        'LegendsNames',{MODEL.esc_names{8}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Post Procesamiento
    % Logaritmo/tendencia vs corr. Anterior
    PostPrLogsComp(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE,'v7', 'PostProcessing'),...
        'Esc',{'v0', MODEL.PostProc},...
        'Esc_add', {'v7', MODEL.PostProc.v7, MODEL.esc_col{7}},...
        'PlotList', pp_list,...
        'LegendsNames',{MODEL.esc_names{8}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Variables en Niveles originales
    PostPrLevels(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v7','PostProcessing'),...
        'Esc_add', {'v7', MODEL.PostProc.v7, MODEL.esc_col{7}},...
        'PlotList', list_lev,...
        'Titles', tit_lev,...
        'LegendsNames',{MODEL.esc_names{8}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Brechas
    PostPrGaps(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v7','PostProcessing'),...
        'Esc_add', {'v7', MODEL.PostProc.v7, MODEL.esc_col{7}},...
        'PlotList', list_gaps,...
        'LegendsNames',{MODEL.esc_names{8}, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Descomposición de choques para variables seleccionadas
    list = {'L_GDP_RW_GAP', 'DLA_CPI_RW', 'RS_RW', 'D4L_CPI_NOSUBY','L_GDP_GAP','DLA_CPIXFE', 'DLA_S', 'D4L_MB', 'RS',...
            'D4L_CPI', 'L_Z_GAP', 'D4L_VEL', 'RR', 'D4L_S', 'RMC', 'MCI'};
    
    % Long
    plot_shd_dsc(MODEL, MODEL.Esc.v7.pred, MODEL.Esc.v7.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v7','Shock_dec\long'),...
        'Rng', {},...
        'Variables',list)

    %short
    plot_shd_dsc(MODEL, MODEL.Esc.v7.pred, MODEL.Esc.v7.shd,...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v7','Shock_dec\short'),...
        'Rng', MODEL.DATES.hist_end-20:MODEL.DATES.hist_end+20,...
        'Variables',list)
    
    % Descomposición de chosques para variables seleccionadas
    % (Primeras diferencias)
    plot_diff_shd_dsc(MODEL,...
                  'variables', list,...
                  'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v7', 'Shock_dec', 'diff'),...
                  'Esc_add', {'v7', MODEL.Esc.v7.shd, MODEL.Esc.v7.pred});
    
    % Contribuciones
    contributions(MODEL,...
        'Esc_add', {'v7', MODEL.Esc.v7.pred}, ... 
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v7', 'contributions') ...           
    );
              
    contributions(MODEL,...
        'Esc_add', {'v7', MODEL.Esc.v7.pred}, ... 
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v7', 'diff_contributions'), ... 
        'Difference', true ...
    );               
    
    % Real exchange rate (subplot)
    tcr_subplot(MODEL,...
        'Esc_add', {'v7', MODEL.Esc.v7.pred, MODEL.esc_col{7}},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.esc_names{8}, MODEL.leg_act});
    
    % Money Velocity (subplot)
    vel_subplot(MODEL,...
        'tab_range', tab_range,...
        'Esc_add', {'v7', MODEL.PostProc, MODEL.esc_col{7}},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.esc_names{8}, MODEL.leg_act},...
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v7', 'otras'));
end
disp('Escenario 7: ok');



