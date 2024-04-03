%%
clear all;
addpath('src')

%% Datos
% Calculo de variables a partir de base de datos inicial
% makedata;
make_obs;
% Calculo de variables no observables a partir de estructura del modelo
% kalmanfilter;

%% parámetros para filtro de kalman
% Período inicial de filtrado
sdate = qq(2001,4);
% Período final de filtrado (fin de historia)
edate = qq(2021,3);

% carga de modelo con p.filter = true (.model)
[m,~,~] = read_QPM(true);
% proceso de filtrado (ModelingLegacy\@model)
[m_kf, g] = filter(m,obs,sdate:edate);
h = g.mean;

%% Tamaño de anclajes
d_rs = max(abs(h.RS{qq(2005,1):end}.diff(-1)));
d_mb = -1*max(abs(h.D4L_MB{qq(2005,1):end}.diff(-1)));

%% Carga de modelo sin filtrado

[m,p,mss] = read_QPM(false);

%% libre
fcstrng = h.L_GDP_GAP.range(end)+1:qq(2030,4);

sl = simulate(m, h, fcstrng, 'anticipate', false, 'DbOverlay=', true);

%% Escenario alterno 

% Base de datos con choques de escenario libre
shocks= sl*get(m,'exList');
db1 = dboverlay(h,shocks);

% Asignación de anclaje
% db1.RS(qq(2021,4)) = 5;
% db1.SHK_RS(qq(2021,4)) = 2.813;%0.795; %
db1.RS(qq(2021,4)) = sl.DLA_CPIXFE(qq(2021,4)) - 2;


% Creación de plan de simulación
planE1 = plan(m,fcstrng);
% Variables a exogenizar y endogenizar
% Dado que estamos dando un choque, no es necesario "endogenizar ninguna
% variable"
planE1 = endogenize(planE1,{'SHK_RS'},qq(2021,4)); 
planE1 = exogenize(planE1,{'DLA_CPIXFE'},qq(2021,4));


% Simulación
s1 = simulate(m,db1,fcstrng,'plan',planE1,'anticipate',false, 'DbOverlay=', true);

%%
hist = databank.fromCSV('output\KalmanHist_long.csv');
bm_ok = databank.fromCSV('data\data_bm.csv','Select=','D4L_MB')
%% GRÁFICAS
plotrng = qq(2015,1):qq(2030,4);


vlist = get(m,'xlist');
vlist{end+1} = 'SHK_D4L_MB';
vlist{end+1} = 'SHK_RS';
% vlist = {'L_GDP_RW_GAP','D4L_CPI_RW','RS_RW',...
%             'D4L_CPI_NOSUBY','D4L_GDP', 'D4L_CPIXFE',...
%             'D4L_S','D4L_MB','RS',...
%             'D4L_CPI', 'D4L_Z','D4L_VEL', 'RR','D4_GDP_SM'};


for i = 1:length(vlist)

    figure('Position', [1 42.0181818181818 1675.63636363636 825.6]);
    n = plot([sl.(vlist{i}){plotrng}, s1.(vlist{i}){plotrng}, ...
              hist.(vlist{i}){plotrng}], ...
              'linewidth',1.5);
          
    set(n(end), 'color','k', 'linewidth',2);
    

    hline(real(mss.(vlist{i})))%,'linestyle','--','linewidth',1.5);
    vline(qq(2021,3));
    title(vlist{i},'fontsize',15,'interpreter','none');
%     set(gca, 'ygrid','on');
    legend({'libre', 'Prueba', 'Observado'}, 'fontsize',13);
    
            saveas(gcf, ...
            fullfile(...
            'graficas\gprueba',...
            strcat(vlist{i},'.png'))...
            );
end
close all
return