%% Ajustes para escenarios
%{
Las variables que pueden verse afectadas en el contexto actual son el
precimiento externo y el IPEI
%}
%% Base de datos inicial y trimestres de ajuste
MODEL.Esc.v0.dbi = MODEL.F;
MODEL.Esc.v0.dates = MODEL.DATES.pred_start:MODEL.DATES.pred_start+3;

%% Brecha del producto externo
MODEL.Esc.v0.dbi.L_GDP_RW_GAP(MODEL.Esc.v0.dates) = [-0.0321, -0.0521, -0.00246, 0.01];
% Anclando la variación cambiaria directamente
% MODEL.Esc.v0.dbi.D4L_S_EP(MODEL.DATES.pred_start:MODEL.DATES.pred_start+7) = [0, 0, 0, 0, 0, 0, 0, 0];

%% Tasa de variación del IPEI
MODEL.Esc.v0.dbi.D4L_IPEI(MODEL.Esc.v0.dates) = [2.8863, 2.9354, 2.85735, 2.7952];
% Anclando valores para la tasa de variación directamente
% MODEL.Esc.v0.dbi.D4S4L_GDP_EP(MODEL.Esc.v0.dates) = [0, 0, 0, 0, 0, 0, 0, 0];

%% Listado de variables a usar en el Plan de Simulación
% Shocks (No cambian)
ListShocksEP = {'SHK_L_GDP_RW_GAP','SHK_D4L_IPEI'};%SHK_D4L_CPI_NOSUBY
% Variables (cambian de acuerdo a lo que se "ancle")
ListVarEP    = {'L_GDP_RW_GAP','D4L_IPEI'};%SHK_D4L_CPI_NOSUBY

%% Creación del Plan de Simulación
% Plan de simulación para escenario Libre (Base)
MODEL.Esc.v0.planSim = plan(MODEL.MF, MODEL.Esc.v0.dates);
% Endogenización de shocks de ajuste
MODEL.Esc.v0.planSim = endogenize(MODEL.Esc.v0.planSim,ListShocksEP,MODEL.Esc.v0.dates); 
% Exógenización de ajustes para variables EP
MODEL.Esc.v0.planSim = exogenize(MODEL.Esc.v0.planSim,ListVarEP,MODEL.Esc.v0.dates);


