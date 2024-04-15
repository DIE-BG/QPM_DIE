%{
Transformación de las variables desde su estado inicial (variables en
frecuencia mensual o trimestral) hasta la trimestralización de variables.

Las principales fuentes son:
    1) FRED (Variables externas)
    2) Banguat (Variables Internas)

Departamento de Investigaciones Económicas - 2024.
MJGM/JGOR

%}
%%
% Creación de folder para datos
if ~isfolder(fullfile('data', 'corrimientos', MODEL.CORR_DATE, MODEL.CORR_VER))
    mkdir(fullfile('data', 'corrimientos', MODEL.CORR_DATE, MODEL.CORR_VER))
end  


%% VARIABLES TRIMESTRALES
% Producto externo y domestico (nominal y real)
% Quarterly
q = databank.fromCSV(fullfile('data', 'raw', MODEL.CORR_DATE, 'quarterly.csv'));

% VARIABLES MENSUALES
% a, a_prom, exp_indx, imp_indx, i_star, cpi_sub, s, bm, cpi
m = databank.fromCSV(fullfile('data', 'raw', MODEL.CORR_DATE, 'monthly.csv'));

%% Construcción de CPI_RW y REM_GDP
% CPI_RW
m.CPI_RW = m.A_prom*m.ind_prec_expus + (1- m.A_prom)*m.ind_prec_impus;

% REM_GDP
q.REM = m.REM.convert('Q', 'method=', @sum);
q.S = m.S.convert('Q', 'method=', 'last');
q.REM_GDP = ((q.REM*q.S)/q.NGDP)*100;

%% Trimestralización
% Se tienen variables stock y de flujo por lo que el proceso de
% trimestralización es disinta para cada una.
% variables stock (trimestralización última del trimestre)
stock = {'CPI', 'CPI_RW', 'CPIXFE', 'RS', 'RS_RW'};

names = dbnames(m);

for i = 1:length(names)
    if ~isempty(strmatch(names{i},stock,'exact'))
        q.(names{i}) = m.(names{i}).convert('Q', 'method', 'last');
    end
end

flow = {'MB'};

for i = 1:length(names)
    if ~isempty(strmatch(names{i},flow,'exact'))
        q.(names{i}) = m.(names{i}).convert('Q', 'method', @mean);
    end
end

% Estructura Preproc (monthly, quarterly, obs)
MODEL.PreProc.monthy = m;
MODEL.PreProc.quarterly = q;
list_mobs = {'GDP', 'CPI', 'CPIXFE', 'S', 'RS', 'GDP_RW', 'CPI_RW', 'RS_RW', 'REM_GDP', 'MB'};
temp = q*list_mobs;

%% Construcción de variables observables para QPM

% Variables en logaritmo
exceptions = {'RS','RS_RW'};

list = dbnames(temp);

for i = 1:length(list)
    if isempty(strmatch(list{i},exceptions,'exact'))
        temp.(['L_' list{i}]) = 100*log(temp.(list{i}));
    end
end

% Brechas HP
% For GDP_RW the observable is L_GDP_RW_GAP 

list = {'L_GDP_RW'};

for i = 1:length(list)
    [temp.([list{i} '_BAR']), temp.([list{i} '_GAP'])] = hpf(temp.(list{i}));
end

% Estructura para observables
MODEL.PreProc.obs = struct();
% Lista de observables (debe coincidir con .model)
obs_list = {'L_GDP', 'L_S','L_CPI','RS','L_CPIXFE','L_CPI_RW','RS_RW', 'L_GDP_RW_GAP', 'L_GDP_RW','REM_GDP', 'L_MB'}; 

for i = 1:length(obs_list)
   MODEL.PreProc.obs.(strcat('OBS_',obs_list{i})) = temp.(obs_list{i}); 
   MODEL.PreProc.obs.(strcat('OBS_',obs_list{i})).UserData.endhist = dat2char(MODEL.DATES.hist_end); 
end

clear temp;

%% Almacenamiento de csv con datos que entran al proceso de filtrado.

databank.toCSV(MODEL.PreProc.obs, MODEL.data_file_name, Inf, 'Decimals=', 5, 'UserDataFields=', 'endhist');