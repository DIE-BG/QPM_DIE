%% Descomposición de choques
%{
    Script genera la descomposición de choques para las variables
    seleccionadas (list), en dos diferentes rangos. 
%}
% Generación de datos
MODEL = Sim.sim.shd_dsc(MODEL);
% Lista de variables a graficar
sh_list = {'L_GDP_RW_GAP', 'DLA_CPI_RW', 'RS_RW', 'D4L_CPI_NOSUBY','L_GDP_GAP','DLA_CPIXFE', 'DLA_S', 'D4L_MB', 'RS',...
            'D4L_CPI', 'L_Z_GAP', 'D4L_VEL', 'RR', 'D4L_S', 'RMC', 'MCI'};

%% Gráficas
% Rango Completo
temp_path_long = fullfile('plots',MODEL.CORR_DATE,MODEL.CORR_VER,'Shock_dec\long');
Sim.scripts.plot_shd_dsc(MODEL, 'SavePath', temp_path_long,...
                             'Rng', {},...
                             'Variables',sh_list);
                         
% Rango corto
temp_path = fullfile('plots',MODEL.CORR_DATE,MODEL.CORR_VER,'Shock_dec\short');
Sim.scripts.plot_shd_dsc(MODEL, 'SavePath', temp_path,...
                             'Rng', MODEL.DATES.hist_end-20:MODEL.DATES.hist_end+20,...
                             'Variables',sh_list); 
                         
temp_path_long = fullfile('plots', MODEL.CORR_DATE, 'v7','Shock_dec');
plot_shd_dsc_comparative(MODEL, MODEL.F_pred, MODEL.Esc.v7.shd, 'SavePath', temp_path_long,...
                             'Rng', MODEL.DATES.hist_end-20:MODEL.DATES.hist_end+20,...
                             'Variables',sh_list,...
                             'Esc', {MODEL_ANT.shd_dsc},...
                             'sub_title', {MODEL.leg_act, 'Agosto 2024'});
