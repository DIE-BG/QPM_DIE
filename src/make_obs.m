%{
DIE - MJGM 01/24

Cálculo de variables observables necesarias para QPM.
A partir de este set de observables se genera la historia del resto de
variables. 

Debe correrse para obtener el resto de variables del modelo a partir de su
estructura (incluyendo tasas de variación anualizadas, interanuales,
brechas, tendencias y otras que no sean observables). 

Genera las variables endógenas y cualquier transformación que esté en
model.model.

La base de dato inicial (data.csv) contiene:
GDP         PIB Trimestral
CPI         IPC Total Trimestral
CPIXFE      Indice Subyacente Optima MSE
S           Tipo de Cambio nominal
RS          Tasa de interés líder
GDP_RW      PIB Estados Unidos
CPI_RW      IPEI
RS_RW       Tasa de fondos federales
REM_GDP     Remesas como % PIB nominal


Variables observables: 

!measurement_variables
OBS_L_GDP,          log(GDP)*100
OBS_L_S,            log(S)*100
OBS_L_CPI,          log(CPI)*100
OBS_RS,             RS
OBS_L_CPIXFE,       log(CPIXFE)*100
OBS_L_CPI_RW,       log(CPI_RW)*100
OBS_RS_RW           RS_RW
OBS_L_GDP_RW_GAP    Brecha filtro hp
OBS_RR_RW_BAR       Tendencia filtro hp
OBS_REM_GDP         REM_GDP

Otras observables que no se le dan al modelo (se recuperan con la
estructura)
(Podrían removerse de la lista en model.model??)

OBS_L_GDP_BAR       empty for filtration    
OBS_L_Z_BAR         empty for filtration
OBS_RR_BAR          empty for filtration


Nota: Debiera agregarse D4L_BM a listado de observables.

%}
%% Load quarterly data
% Base de datos inicial
d = databank.fromCSV('data/data.csv');
% Estructura para cálculos intermedios
dd = struct();
%% Make log of variables
exceptions = {'RS','RS_RW'};

list = dbnames(d);

for i = 1:length(list)
    if isempty(strmatch(list{i},exceptions,'exact'))
        dd.(['L_' list{i}]) = 100*log(d.(list{i}));
    end
end

dd.RS = d.RS;
dd.RS_RW = d.RS_RW;
dd.REM_GDP = d.REM_GDP;

%% Real Interest Rates
% Foreign real interest rate (history)
dd.RR_RW = d.RS_RW - (dd.L_CPI_RW - dd.L_CPI_RW{-4});

%% Trends and Gaps - Hodrick-Prescott filter 
% For GDP_RW and RR_RW the observables for the QPM are L_GDP_RW_GAP and
% RR_RW_BAR

list = {'RR_RW','L_GDP_RW'};

for i = 1:length(list)
    [dd.([list{i} '_BAR']), dd.([list{i} '_GAP'])] = hpf(dd.(list{i}));
end

%% select just observables
% Estructura para observables
obs = struct();
% Lista de observables (debe coincidir con .model)
obs_list = {'L_GDP', 'L_S','L_CPI','RS','L_CPIXFE','L_CPI_RW','RS_RW','L_GDP_RW_GAP','REM_GDP','L_MB'}; 

for i = 1:length(obs_list)
   obs.(strcat('OBS_',obs_list{i})) = dd.(obs_list{i}); 
end

% Otras observables
% obs.OBS_L_GDP_BAR    = tseries();
% obs.OBS_L_Z_BAR      = tseries();
% obs.OBS_RR_BAR       = tseries();

%% Almacenamiento de base de datos
databank.toCSV(obs,'output\observables.csv', Inf);

