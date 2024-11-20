% Datos Full data QPM
disp('-----> Carga de datos QPM');
qpm_act = databank.fromCSV(fullfile('data', 'fulldata', MODEL.CORR_DATE, 'fulldata.csv'));

%% MENSUALIZACIÓN
  disp('-----> Mensualizacion de Inflacion Subyacente Óptima MSE');
  mens = struct;
  mens.v_cpi_sub = tseries(mmtoday:mm(2030,12), 0);
  mens.v_cpi_sub = convert(qpm_act.D4L_CPIXFE,12,'method=','linear','position=','end');
  mens.v_cpi_sub = resize(mens.v_cpi_sub, mmtoday:mm(2030,12));

  disp('-----> Mensualizacion de Inflacion Total del QPM');
  mens_2 = struct;
  mens_2.v_cpi = tseries(mmtoday:mm(2030,12), 0);
  mens_2.v_cpi = convert(qpm_act.D4L_CPI,12,'method=','linear','position=','end');
  mens_2.v_cpi = resize(mens_2.v_cpi, mmtoday:mm(2030,12));

%% Almacenamiento de Inlación Subaycente e Inflación Total Mensual del QPM.
databank.toCSV(mens,fullfile('data', 'fulldata', MODEL.CORR_DATE, 'Inf_sub_mensuales.csv'), Inf);
databank.toCSV(mens_2,fullfile('data', 'fulldata', MODEL.CORR_DATE, 'Inf_tot_mensuales.csv'), Inf);
disp('Almacenamiento terminado');

%% Datos Promedio (QPM y MME)
% disp('-----> Carga de datos promedio');
% prom_act = databank.fromCSV(fullfile('data', 'fulldata', MODEL.CORR_DATE, 'datos_esc_base_nov24.csv'));
% 
% %% MENSUALIZACIÓN PROMEDIO
%   disp('-----> Mensualizacion de Inflacion Total Promedio');
%   mens_3 = struct;
%   mens_3.v_cpi_prom = tseries(mmtoday:mm(2030,12), 0);
%   mens_3.v_cpi_prom = convert(prom_act.inf_prom,12,'method=','linear','position=','end');
%   mens_3.v_cpi_prom = resize(mens_3.v_cpi_prom, mmtoday:mm(2030,12));
% 
% %% Almacenamiento de Inflación Total Promedio Mensual.
% databank.toCSV(mens_3,fullfile('data', 'fulldata', MODEL.CORR_DATE, 'Inf_promedio_mensuales.csv'), Inf);
% disp('Almacenamiento terminado');


