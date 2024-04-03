%%%%%%%%%
%%% FILTRATION
%%%%%%%%%


%% Data sample (to extract non observables)
sdate = qq(2001,4);
edate = qq(2023,3);



%% Load data
d = dbload('output\history.csv');

% db con observables para filtrado
% observables base de datos 
%(todas excepto GDP_RW, la observable es la brecha)
% Estas deben coindidir con las de model.model
dd.OBS_L_GDP        = d.L_GDP;
dd.OBS_L_CPI        = d.L_CPI;
dd.OBS_L_CPIXFE     = d.L_CPIXFE;
dd.OBS_L_S          = d.L_S;
dd.OBS_RS           = d.RS;
dd.OBS_L_CPI_RW     = d.L_CPI_RW;
dd.OBS_RS_RW        = d.RS_RW;
dd.OBS_REM_GDP      = d.REM_GDP;


% Otras observables calculadas en make data (filtro HP) y vacías
% Autorregresivo en QPM, mejor darle la observada
dd.OBS_L_GDP_RW_GAP = d.L_GDP_RW_GAP;
% Autorregresivo en QPM
dd.OBS_RR_RW_BAR    = d.RR_RW_BAR;

% Se generan a partir del proceso de filtrado
dd.OBS_L_GDP_BAR    = tseries();
dd.OBS_L_Z_BAR      = tseries();
dd.OBS_RR_BAR       = tseries();

%% Reads the model (filter=true for conditional in model.model, 
% must be a parameter within the structure assigned to the model object)

[m,p,mss] = read_QPM(true);

% Filtration (filter function from IRIS Toolbox, not matlab)
%{
m_kf: filtered model
g: Output database or tseries object with the DFM observables.
%}
[m_kf, g] = filter(m,dd,sdate:edate);
h = g.mean;
dk = dbextend(d,h);

%% Save the database
databank.toCSV(dk,'output\kalm_his.csv', Inf);
