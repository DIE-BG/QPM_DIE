%{
Creación del Escenario en el que se modifica el valor del coeficiente
autorregresivo asociado a las remesas como porcentaje del PIB.
%}

% parámetros del QPM
setparam;
s.rho_REM_GDP = 0.90;

% Modelo
M = model('QPM.model', 'assign', s);
M = sstate(M,'growth=',true,'MaxFunEvals',1000,'display=','off');
M = solve(M,'error=',true);
    
MODEL.Esc.v6.MODEL = M;

% Filtro de kalman
[MODEL.Esc.v6.MF, MODEL.Esc.v6.F] = filter(MODEL.Esc.v6.MODEL, MODEL.PreProc.obs,... 
                                MODEL.DATES.hist_start:MODEL.DATES.hist_end, ... 
                                'meanOnly=',true);
                            
% Pronosticos
fcstrng = MODEL.DATES.pred_start:MODEL.DATES.pred_end;
MODEL.Esc.v6.pred = simulate(MODEL.Esc.v6.MF,... Modelo Filtrado
                        MODEL.Esc.v6.F,...Base de datos inicial filtrada
                        fcstrng,... Rango de Simulación
                        'anticipate', false,...
                        'DbOverlay=', true);

% Descopmposición de shocks
MODEL.Esc.v6.shd = simulate(MODEL.Esc.v6.MF,...
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

%% Cambio de etiqueta para el modelo con media movil
leg_act = 'Octubre 2024'; 
MODEL.Esc.v6.name = 'Octubre 2024 rho_rem_gdp';

%% Graficas
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
        'Esc_add', {'v6', MODEL.Esc.v6.pred, MODEL.esc_col{6}},...
        'PlotList', get(MODEL.MF, 'xlist'),...
        'LegendsNames',{MODEL.Esc.v6.name, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    
    % Post Procesamiento
    % Logaritmo/tendencia vs corr. Anterior
    PostPrLogsComp(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE,'v6', 'PostProcessing'),...
        'Esc',{'v0', MODEL.PostProc},...
        'Esc_add', {'v6', MODEL.PostProc.v6, MODEL.esc_col{6}},...
        'PlotList', pp_list,...
        'LegendsNames',{MODEL.Esc.v6.name, MODEL.leg_act},...
        'TabRange', tab_range...
        );
    % Variables en Niveles originales
    PostPrLevels(MODEL,...
        'StartDate',{MODEL.DATES.hist_start, MODEL.DATES.hist_end - 20},...
        'EndDatePlot', {MODEL.DATES.pred_end, MODEL.DATES.hist_end + 20},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v6','PostProcessing'),...
        'Esc_add', {'v6', MODEL.PostProc.v6, MODEL.esc_col{6}},...
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
        'Esc_add', {'v6', MODEL.PostProc.v6, MODEL.esc_col{6}},...
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
        'Esc_add', {'v6', MODEL.Esc.v6.pred, MODEL.esc_col{6}},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.Esc.v6.name, MODEL.leg_act});
    
    % Money Velocity (subplot)
    vel_subplot(MODEL,...
        'tab_range', tab_range,...
        'Esc_add', {'v6', MODEL.PostProc, MODEL.esc_col{6}},...
        'tab_range', tab_range,...
        'LegendsNames',{MODEL.Esc.v6.name, MODEL.leg_act},...
        'SavePath', fullfile('plots', MODEL.CORR_DATE, 'v6', 'otras'));
end
disp('Escenario 6: ok');