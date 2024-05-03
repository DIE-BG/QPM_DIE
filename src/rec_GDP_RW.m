function [L_GDP_RW, L_GDP_RW_BAR, D4L_GDP_RW, DLA_GDP_RW, D4_GDP_RW_SM] = rec_GDP_RW(q, F_pred, d, varargin)

%{
Se realiza la reconstrucción de GDP_RW partiendo de pornósticos . 

## Syntax ##

    MODEL = PpostProcessing(MODEL, varargin)

## Input Arguments ##

__`MODEL`__ [ struct ] -
Debe contener al menos la estructura con los resultados del proceso de
simulación MODEL.F_pred.

- DIE
- Marzo 2024
- MJGM/JGOR
%}

p = inputParser;
    addParameter(p, 'list', {});
    addParameter(p, 'list_niv', {});
parse(p, varargin{:});
params = p.Results; 

MODEL.PreProc.quarterly = q;
MODEL.F_pred = F_pred;
MODEL.DATES = d;

%% Modelo ARIMA(1,1,0)
md = arima(1,1,0);
estmd = estimate(md, MODEL.PreProc.quarterly.L_GDP_RW_BAR.Data);

fcstrng = length(MODEL.PreProc.quarterly.L_GDP_RW.End+1:MODEL.DATES.pred_end);

[Y, ~] = forecast(estmd,  fcstrng,...
                    MODEL.PreProc.quarterly.L_GDP_RW_BAR.Data);

MODEL.PreProc.quarterly.L_GDP_RW_BAR(MODEL.PreProc.quarterly.L_GDP_RW.End+1:MODEL.DATES.pred_end) = Y;
MODEL.F_pred.L_GDP_RW_BAR = MODEL.PreProc.quarterly.L_GDP_RW_BAR;

% GDP_RW
MODEL.F_pred.L_GDP_RW = MODEL.F_pred.L_GDP_RW_GAP + MODEL.F_pred.L_GDP_RW_BAR;

% D4L_GDP_RW
MODEL.F_pred.D4L_GDP_RW = MODEL.F_pred.L_GDP_RW.diff(-4);
% DLA_GDP_RW
MODEL.F_pred.DLA_GDP_RW = MODEL.F_pred.L_GDP_RW.diff(-1)*4;

% suma movil
temp_movsum = movsum(MODEL.F_pred.L_GDP_RW.Data, 4);
temp_movsum = temp_movsum(3:end-1);
MODEL.F_pred.GDP_RW_SM = MODEL.F_pred.L_GDP_RW;
MODEL.F_pred.GDP_RW_SM.Data = temp_movsum;

% Tasa de variación de la suma movil
MODEL.F_pred.D4_GDP_RW_SM = MODEL.F_pred.GDP_RW_SM.diff(-4)/4; 

% Nombres de las variables
MODEL.F_pred.L_GDP_RW.comment = 'Producto Interno Bruto Real EEUU (Logaritmo)';
MODEL.F_pred.L_GDP_RW_BAR.comment= 'Tendencia del Producto Interno Bruto Real EEUU';
MODEL.F_pred.D4L_GDP_RW.comment = 'Tasa de Variación Internual del PIB de EEUU';
MODEL.F_pred.DLA_GDP_RW.comment = 'Tasa de Variación Intertrimestral anualizada del PIB de EEUU';
MODEL.F_pred.D4_GDP_RW_SM.comment = 'Tasa de Variación Interanual de la suma de 4 Trimestres del PIB de EEUU';


L_GDP_RW = MODEL.F_pred.L_GDP_RW;
L_GDP_RW_BAR = MODEL.F_pred.L_GDP_RW_BAR;
D4L_GDP_RW = MODEL.F_pred.D4L_GDP_RW;
DLA_GDP_RW = MODEL.F_pred.DLA_GDP_RW;
D4_GDP_RW_SM = MODEL.F_pred.D4_GDP_RW_SM;

end
