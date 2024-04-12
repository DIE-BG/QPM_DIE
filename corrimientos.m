%%
%{
QPM
%}

clear all;
addpath('src')

%% Carga de configuraciones generales del corrimiento base (v0)
StartConfig;

%% Preprocessing
Preprocessing;
disp('Preprocesamiento: ok');
%% Observables
% Calculo de variables a partir de base de datos inicial
% makedata;
make_obs;
disp('Makeobs: ok');
%% parámetros para filtro de kalman
% Período inicial de filtrado
sdate = qq(2001,4);
% Período final de filtrado (fin de historia)
edate = qq(2023,3);

% carga de modelo con p.filter = true (.model)
[m,~,~] = read_QPM(true);
% proceso de filtrado (ModelingLegacy\@model)
[m_kf, g] = filter(m,obs,sdate:edate);
h = g.mean;

%% Carga de modelo sin filtrado

[m,p,mss] = read_QPM(false);

%% libre
fcstrng = h.L_GDP_GAP.range(end)+1:qq(2030,4);

sl = simulate(m, h, fcstrng, 'anticipate', false, 'DbOverlay=', true);


%%


subplot(2,2,1)
plot([sl.L_Z_GAP]); zeroline;

subplot(2,2,2);
plot([sl.D4L_CPI, sl.D4L_CPIXFE]); hline(4);

subplot(2,2,3)
plot(exp(sl.L_S/100)); hline(mean(exp(sl.L_S/100)));

subplot(2,2,4)
plot([sl.DLA_CPI_RW, sl.D4L_CPI_RW]); zeroline;


%%
figure;
subplot(2,2,1)
plot(sl.RMC);
title('RMC');
zeroline;

subplot(2,2,2)
plot([sl.L_GDP_GAP, sl.L_Z_GAP]);
title('RMC COMPONENTS');
legend({'y_gap', 'z_gap'}, 'Interpreter','none');
zeroline;

subplot(2,2,3)
plot(sl.MCI); 
title('MCI');
zeroline;

subplot(2,2,4)
plot([sl.RR_GAP, sl.L_Z_GAP, sl.D4L_MB]);
title('MCI COMPONENTS');
legend({'RR_GAP', 'z_gap', 'MB'}, 'Interpreter','none');
zeroline;