%{
Función para asignación de parametros de modelo QPM
Original: curso MPAF IMF
Primera versión: Ajustes por Juan Pablo Medina

Esta versión: DIE-MJGM 01/2024

%}
function [m,p,mss] = read_QPM(filter)

%% Filtration on/off
% filter = true - Kalman filter ON
% filter = false - Kalman filter OFF
p.filter = filter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Aggregate demand equation (the IS curve)
% L_GDP_GAP = b1*L_GDP_GAP{-1} - b2*MCI + b3*L_GDP_RW_GAP + ...
%     b5*(REM_GDP-ss_REM_GDP) + SHK_L_GDP_GAP;

% Real monetary conditions index (mci)
% MCI = b4*RR_GAP + (1-b4)*(-L_Z_GAP) - b6*(D4L_MB{-1}-ss_D4L_MB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%output persistence;
%b1 varies between 0.1 (extremely flexible) and 0.95(extremely persistent)
p.b1 = 0.3526;%0.4310;

%policy passthrough (the impact of monetary policy on real economy); 
%b2 varies between 0.1 (low impact) to 0.5 (strong impact)
p.b2 = 0.1422;%0.1572; 

%the impact of external demand on domestic output; 
%b3 varies between 0.1 and 0.7
p.b3 = 0.558;%0.5905; 

%the weight of the real interest rate and real exchange rate gaps in MCI;
%b4 varies from 0.3 to 0.8
p.b4 = 0.8061;%0.8013;

%Direct effect of remmittances on output gap (after controlling by foreign GDP) 
p.b5 = 0.1183;%0.0702;

%Effect of money growth in MCI
p.b6 = 1; %0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Aggregate supply equation (the Phillips curve)
% DLA_CPIXFE =  a1*DLA_CPIXFE{-1} + (1-a1)*(DLA_CPIXFE{1}) + a2*RMC + ...
%               SHK_DLA_CPIXFE;

% Real marginal cost (rmc)
% RMC = a3*L_GDP_GAP + (1-a3)*(L_ZIMP_GAP);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inflation persistence; 
% a1 varies between 0.4 (implying low persistence) to 0.9 
% (implying high persistence)
p.a1 = 0.3707;%0.4053; 

% policy passthrough (the impact of rmc on inflation); 
% a2 varies between 0.1 (a flat Phillips curve and a high sacrifice ratio) 
% to 0.5 (a steep Phillips curve and a low sacrifice ratio)
p.a2 = 0.05;%0.017;%0.0176;

% the ratio of imported goods in firms' marginal costs (1-a3); 
% a3 varies between 0.9 for a closed economy to 0.5 for an open economy
p.a3 = 0.5;%0.8193;%0.7396;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Uncovered Interest Rate Parity (UIP)
% B. dirty float/managed float/peg
% L_S = h2*(L_S{-1} + DLA_S_TAR/4) + (1-h2)*((1-e1)*L_S{+1} +e1*(L_S{-1} +...
%       2/4*(D4L_CPI_TAR - ss_DLA_CPI_RW + DLA_Z_BAR)) + ...
%       (- RS + h3*(D4L_MB - ss_D4L_MB) + RS_RW + PREM)/4) + SHK_L_S;

% UIP when the exchange rate is managed to meet the inflation objective
% Set h2 different from zero. See discussion related to Section B and
% parameter h2. On the top of that set f2 and f2 different from zero.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting e1 equal to 0 reduces the equation to the simple UIP
p.e1 = 0.3988;%0.4948;

% Setting h2 equal to 0 implies float based on the UIP, h2 equal to 1 means
% managed FX following target appreciation.
p.h2 = 0.822;%0.7256;

% D. Exchange Rate Target
% DLA_S_TAR = f1*DLA_S_TAR{-1} + (1-f1)*(D4L_CPI_TAR - ss_DLA_CPI_RW + ...
%             DLA_Z_BAR + f2*(D4L_CPI-D4L_CPI_TAR) + f3*L_GDP_GAP) + SHK_DLA_S_TAR;
% L_S_TAR = L_S_TAR{-1} + DLA_S_TAR/4;
p.f1 = 0.845;%0.7437;
p.f2 = 0;
p.f3 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Monetary policy reaction function 

% The rule modified for imperfect control over the domestic money market
%(use if the central bank stabilizes the exchange rate by FOREX interventions)
% RS = h1*(4*(L_S{+1} - L_S) + RS_RW + PREM) + (1 - h1)*(g1*RS{-1} + ...
%             (1 - g1)*(RSNEUTRAL + g2*(D4L_CPI{+4} - D4L_CPI_TAR{+4}) + ...
%             g3*L_GDP_GAP + g4*(DLA_S{-1}-DLA_S_TAR{-1}))) + SHK_RS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% policy persistence; 
% g1 varies from 0 (no persistence) to 0.8 ("wait and see" policy)
p.g1 = 0.9578;%0.9832; 

% policy reactiveness: the weight put on inflation by the policy maker); 
% g2 has no upper limit but must be always higher then 0 (the Taylor principle)
p.g2 = 0.5; 

% policy reactiveness: the weight put on the output gap by the policy maker); 
% g3 has no upper limit but must be always higher then 0
p.g3 = 0.4932;%0.4520;

% degree to which the central bank does not control domestic money market
p.h1 = 0.0;

p.h3 = 0.25;
% reaction to exchange rate devaluation with respect its target
p.g4 = 0; %1.7334;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Money demand equation for monetary base or other monetary aggregate'
%D4L_MB = D4L_CPI + j1*D4L_GDP - j2*(RS-RS{-4}) - D4L_VEL;
p.j1 = 1.0;
p.j2 = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Speed of convergence of selected variables to their trend values.
% Used for risk premium, trends, foreign variables and commodity prices, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% persistent shock to risk premium
% SHKN_PREM = rho_SHKN_PREM*SHKN_PREM{-1} + SHK_PREM;
p.rho_SHKN_PREM = 0.5;

% persistence in convergence of trend variables to their steady-state levels
% applies for:   DLA_GDP_BAR, DLA_Z_BAR, RR_BAR and RR_RW_BAR

p.rho_DLA_Z_BAR   = 0.8368;%0.7523;
p.rho_DLA_GDP_BAR = 0.8296;%0.7731;
p.rho_RR_BAR      = 0.8452;%0.7662;
p.rho_RR_RW_BAR   = 0.852;%0.8614;

% persistence in foreign GDP 
% L_GDP_RW_GAP = h2*L_GDP_RW_GAP{-1} + SHK_L_GDP_RW_GAP;
p.rho_L_GDP_RW_GAP = 0.8;%0.7806;
% adding persistence of remittance shock 
p.rho_REM_GDP = 0.9894;%0.9883;
% adding persistence shock to D4L_CPI_NOSUBY 
p.rho_D4L_CPI_NOSUBY = 0.7893;%0.7823;
%adding persistence shock to money velocity
p.rho_D4L_VEL = 0.8528;%0.8513;
p.rho_PM_D4L_MB = 0.7;

% persistence in foreign interest rates (and inflation);
%RS_RW = rho_RS_RW*RS_RW{-1} + (1-rho_RS_RW)*(RR_BAR + DLA_CPI_RW) + SHK_RS_RW;
p.rho_RS_RW      = 0.9876;%0.8742;
p.rho_DLA_CPI_RW = 0.7;%0.4724;%0.5666;

% Speed of inflation target adjustment to the medium-term target (higher values mean slower adjustment)
% D4L_CPI_TAR = f1*D4L_CPI_TAR{-1} + (1-f1)*ss_D4L_CPI_TAR + SHK_D4L_CPI_TAR;
p.rho_D4L_CPI_TAR = 0.8723;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibrated Steady States
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Foreign trend inflation or inflation target
p.ss_DLA_CPI_RW = 2;
% Trend level of domestic real interest rate 
p.ss_RR_BAR = 1.0; %0.5; 
% Trend change in the real ER (negative number = real appreciation)
p.ss_DLA_Z_BAR = -2.0; 
% Potential output growth
p.ss_DLA_GDP_BAR = 3.5; %1
% Trend level of foreign real interest rate
p.ss_RR_RW_BAR = 0.50; %0.5;
% Domestic inflation target
p.ss_D4L_CPI_TAR = 4;  
% Steady state value for Remittance to GDP ratio
p.ss_REM_GDP  = 15.0;
%adding money growth steady state
p.ss_D4L_MB = 10;

p.ss_D4L_CPIXFE = 4;
p.ss_D4L_VEL = -2.37;
p.ss_D4L_GDP = 3.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard Deviation for Shocks in QPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.std_SHK_L_GDP_GAP      = 1; %estimation: 1.2035
p.std_SHK_DLA_GDP_BAR = 0.6; %estimation: 0.8201

p.std_SHK_D4L_CPI_TAR = 0; %estimation: 0.0046

p.std_SHK_DLA_CPIXFE     = 0.9; %estimation: 0.9935
p.std_SHK_D4L_CPI_NOSUBY = 0.9; %estimation: 1.0522

p.std_SHK_D4L_VEL  = 0.5;%2; %estimation: 1.9985
p.std_SHK_D4L_MB   = 4.6711; 

p.std_SHK_L_S      = 1.1784;%3; %estimation: 0.9710
p.std_SHK_DLA_S_TAR   = 0.5;%1; %estimation: 1.9413
p.std_SHK_DLA_Z_BAR   = 0.75; %estimation: 1.6697

p.std_SHK_PREM         = 1;

p.std_SHK_RS         = 0.4572;%1; %estimation: 0.4868
p.std_SHK_RR_BAR      = 0.3; %estimation: 0.1854

p.std_SHK_REM_GDP       = 0.5; %estimation: 0.5344
p.std_SHK_L_GDP_RW_GAP  = 1; %estimation: 1.3091
p.std_SHK_RS_RW         = 0.4713;%1; %estimation: 0.2952
p.std_SHK_DLA_CPI_RW    = 9.7052;%estimation: 9.6428
p.std_SHK_RR_RW_BAR     = 0.25;%0.5; %estimation: 1.5696

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model solving
m = model('model.model','linear=',false,'assign',p);
m = sstate(m,'growth',true,'MaxFunEvals',2000,'display=','off');
mss = get(m,'sstate');

%% Check steady state
[flag,discrep,eqtn] = chksstate(m);

if ~flag
  error('Equation fails to hold in steady state: "%s"\n', eqtn{:});
end

m = solve(m);
end
