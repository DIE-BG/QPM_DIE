%% Ajustes para variables EP
%{
    Las variables de ajuste para presentaciones EP son
    Crecimiento Económico
    - D4S4L_GDP_ADJ: Ajustes
    - D4S4L_GDP_EP: Variable final a presentar
    Depreciación Cambiaria interanual y Tipo de cambio relacionado
    - D4L_S_ADJ: Ajustes
    - D4L_S_EP: Variable final a presentar
    - L_S_EP: Reconstrucción del nivel ajustado
    Tasa de interés líder de política monetaria
    - RS_ADJ: Ajustes.
    - i_EP: Recomendación a presentar.
    Inflación Total Interanual
    - D4L_CPI_NOSUBY: Inflación NO subyacente, para ajustes.
    - D4L_CPI: Inflación total. Suma de la D4L_CPIXFE + D4L_CPI_NOSUBY.

Los ajustes pueden hacerse de dos formas:
    1. Anclando valores para los ajustes:
        - Se exogeniza la variable "_ADJ" 
        - Se endogeniza el shock de la misma variable "SHK_*_ADJ"
          (SHK_D4L_CPI_NOSUBY en el caso de la inflación)
    2. Anclando los valores finales de la variable a presentar
        - Se exógeniza la variable "_EP"
        - Se endogeniza el shock de la variable de ajuste "SHK_*_ADJ"
          (SHK_D4L_CPI_NOSUBY en el caso de la inflación)
          (igual que en el caso 1).
%}
%% Base de datos inicial y trimestres de ajuste
MODEL.Esc.v0.dbi = MODEL.F;
MODEL.Esc.v0.dates = MODEL.DATES.pred_start:MODEL.DATES.pred_start+7;

%% Tasa de interés líder (VARIABLE PARA RECOMENDACIÓN: i_EP)
% Anclado Ajustes:
MODEL.Esc.v0.dbi.RS_ADJ(MODEL.Esc.v0.dates) = [0, 0, 0, 0, 0, 0, 0, 0];
% Anclando valores para la recomendación
% MODEL.Esc.v0.dbi.i_EP(MODEL.Esc.v0.dates) = [5.25, 5.5, 5.75, 6, 6, 6, 6, 6];

%% Inflación total interanual (VARIABLE FINAL PARA EP: D4L_CPI)
% Anclando Ajustes
% MODEL.Esc.v0.dbi.D4L_CPI_NOSUBY(MODEL.DATES.pred_start:MODEL.DATES.pred_start+7) = [0, 0, 0, 0, 0, 0, 0, 0];
% Anclando inflación directamente
% MODEL.Esc.v0.dbi.D4L_CPI(MODEL.DATES.pred_start:MODEL.DATES.pred_start+7) = [0, 0, 0, 0, 0, 0, 0, 0];

%% Tasa de variación del tipo de cambio nominal (VARIABLES FINALES: D4L_S_EP y L_S_EP)
MODEL.Esc.v0.dbi.D4L_S_ADJ(MODEL.Esc.v0.dates) = [0, 0, 0, 0, 0, 0, 0, 0];
% Anclando la variación cambiaria directamente
% MODEL.Esc.v0.dbi.D4L_S_EP(MODEL.DATES.pred_start:MODEL.DATES.pred_start+7) = [0, 0, 0, 0, 0, 0, 0, 0];

%% Tasa de variación de la suma móvil de 4 trimestres del PIB (VARIABLE FINAL: D4S4L_GDP_EP)
MODEL.Esc.v0.dbi.D4S4L_GDP_ADJ(MODEL.Esc.v0.dates) = [0, 0, 0, 0, 0, 0, 0, 0];
% Anclando valores para la tasa de variación directamente
% MODEL.Esc.v0.dbi.D4S4L_GDP_EP(MODEL.Esc.v0.dates) = [0, 0, 0, 0, 0, 0, 0, 0];

%% Listado de variables a usar en el Plan de Simulación
% Shocks (No cambian)
ListShocksEP = {'SHK_D4S4L_GDP_ADJ','SHK_D4L_S_ADJ','SHK_RS_ADJ'};%SHK_D4L_CPI_NOSUBY
% Variables (cambian de acuerdo a lo que se "ancle")
ListVarEP    = {'D4S4L_GDP_ADJ','D4L_S_ADJ','RS_ADJ'};%SHK_D4L_CPI_NOSUBY

%% Creación del Plan de Simulación
% Plan de simulación para escenario Libre (Base)
MODEL.Esc.v0.planSim = plan(MODEL.MF, MODEL.Esc.v0.dates);
% Endogenización de shocks de ajuste
MODEL.Esc.v0.planSim = endogenize(MODEL.Esc.v0.planSim,ListShocksEP,MODEL.Esc.v0.dates); 
% Exógenización de ajustes para variables EP
MODEL.Esc.v0.planSim = exogenize(MODEL.Esc.v0.planSim,ListVarEP,MODEL.Esc.v0.dates);


