%{
DIE - MJGM 01/24

Proceso de filtrado con filtro de Kalman a partir de variables observables.

Se utiliza para recuperar la base de datos histórica que se utiliza 
(o podría utilizarse) como base de condiciones inciales para el proceso
de generación de pronósticos posterior. 

(La base de datos de condiciones iniciales puede provenir de otros 
procedimientos como filtros univariados, fuentes externas, etc.)

Partiendo de data.csv y make_obs.m, se genera un set de 13 variables 
observables (10 para efectos prácticos). 

A partir de estas 13 (10) variables, se genera el resto de variables del 
modelo utilizando la estructura de model.model. Esto incluye brechas, tasas
de variación interanual, intertrimestral anualizadas, brechas, tendencias y
otras.

%}
% clear all;
addpath('src')
%% Datos
% Calculo de variables observables a partir de base de datos inicial
make_obs;

%% parámetros para filtro de kalman
% Período inicial de filtrado
sdate = qq(2001,4);
% Período final de filtrado (fin de historia)
edate = qq(2023,3);

% carga de modelo con p.filter = true (.model)
[m,p,mss] = read_QPM(true);
% proceso de filtrado (ModelingLegacy\@model)
[m_kf, g] = filter(m,obs,sdate:edate);
h = g.mean;

%% Save the database
databank.toCSV(h,'output\KalmanHist.csv', Inf);