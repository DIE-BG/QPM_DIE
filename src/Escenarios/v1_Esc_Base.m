%% Escenario Base para EP
%{
    Se realizan los ajustes para las variables principales
    Crecimiento Económico
    - D4_GDP_SM_ADJ
    - D4_GDP_SM_EP
    Inflación Total
    - D4L_CPI_ADJ
    - D4L_CPI_EP
    Depreciación Cambiaria interanual y Tipo de cambio relacionado
    - D4L_S_ADJ
    - D4L_S_EP
    - L_S_EP
    Tasa de interés líder de política monetaria
    - i_adj
    - i_EP
%}

%% Ajustes Escenario Base para EP

% EJEMPLO!! 
% Fechas en las que se harán ajustes (horizonte de gráficas principales?)
MODEL.DATES.E1_dates = MODEL.DATES.pred_start:qq(2025,4);

% Base de datos inicial
MODEL.Esc.v1.dbi = databank.clip(MODEL.F, MODEL.F_pred.L_CPIXFE.Start, MODEL.F_pred.L_CPIXFE.End);

% Ajustes para "dibujar" trayectorias
% Mantener la tasa en 5%  EN TODA la trayectoria
MODEL.Esc.v1.dbi.i_adj(MODEL.DATES.E1_dates) = 5.00 - MODEL.F_pred.RS((MODEL.DATES.E1_dates));
% Trayectoria de inflación 0.2 puntos por encima EN TODA la trayectoria
MODEL.Esc.v1.dbi.D4L_CPI_ADJ(MODEL.DATES.E1_dates) = 0.2;

% Variables que no es necesario ajustar (por ejemplo)
MODEL.Esc.v1.dbi.D4_GDP_SM_ADJ(MODEL.DATES.E1_dates) = 0;
MODEL.Esc.v1.dbi.D4L_S_ADJ(MODEL.DATES.E1_dates) = 0;

% Lista de variables que ajustan:
listShAdj = {'SHK_i_adj','SHK_D4L_CPI_ADJ','SHK_D4_GDP_SM_ADJ','SHK_D4L_S_ADJ'};
listVarsAdj = {'i_adj','D4L_CPI_ADJ','D4_GDP_SM_ADJ','D4L_S_ADJ'};

% Plan de simulación para una trayectoria de choques definida en los
% trimestres de MODEL.Esc.E1_dates.
MODEL.Esc.v1.planSim = plan(MODEL.MF, MODEL.DATES.pred_start:MODEL.DATES.pred_end);
% Variable a endogenizar (shock propio?? No necesariamente)
MODEL.Esc.v1.planSim = endogenize(MODEL.Esc.v1.planSim,listShAdj,MODEL.DATES.E1_dates); 
% Variable a exogenizar (Anclaje)
MODEL.Esc.v1.planSim = exogenize(MODEL.Esc.v1.planSim,listVarsAdj,MODEL.DATES.E1_dates);

% Simulación.
MODEL.Esc.v1.pred = simulate(MODEL.MF,...
                  MODEL.Esc.v1.dbi,...
                  qq(2022,1):MODEL.DATES.pred_end,...
                  'plan',MODEL.Esc.v1.planSim,...
                  'anticipate',false,...
                  'DbOverlay=', true);

listGraph = {'D4_GDP_SM_EP','D4L_CPI_EP','D4L_S_EP','i_EP'};

%% Post-Procesamiento de variables seleccionadas.
pp_list = {'L_MB', 'L_VEL', 'L_CPI_RW', 'L_CPI_RW_Q','L_Z', 'L_GDP', 'L_GDP_RW'};
list_nivel = {'L_S','L_MB'};

[MODEL.Esc.v1.pred.L_GDP_RW, MODEL.Esc.v1.pred.L_GDP_RW_BAR,...
    MODEL.Esc.v1.pred.D4L_GDP_RW, MODEL.Esc.v1.pred.DLA_GDP_RW,...
    MODEL.Esc.v1.pred.D4_GDP_RW_SM] = rec_GDP_RW(databank.clip(MODEL.PreProc.quarterly, MODEL.DATES.hist_start, qq(2021,4)),...
                                            MODEL.Esc.v1.pred,...
                                            MODEL.DATES);

MODEL = PostProcessing(MODEL,...
    'list',pp_list,...
    'list_niv', list_nivel,...
    'Esc',{'v1', MODEL.Esc.v1.pred});

%% Graficas

    simPlots(MODEL,...
        'StartDate',{qq(2015,4)},...
        'EndDatePlot', {qq(2025,4)},...
        'SavePath', fullfile(cd, 'plots', MODEL.CORR_DATE, 'v1','prediction_compared'),...}
        'Esc_add', {'v1', MODEL.Esc.v1.pred},...
        'PlotList', listGraph,...%get(MODEL.MF, 'xlist'),...
        'LegendsNames',{'Esc. Base', MODEL.leg_act},...
        'TabRange', tab_range...
        );
 