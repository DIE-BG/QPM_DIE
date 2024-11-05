s = struct();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Aggregate demand equation (the IS curve)
% L_GDP_GAP = b1*L_GDP_GAP{-1} - b2*MCI + b3*L_GDP_RW_GAP + ...
%     b5*(REM_GDP-ss_REM_GDP) + SHK_L_GDP_GAP;

% Real monetary conditions index (mci)
% MCI = b4*RR_GAP + (1-b4)*(-L_Z_GAP) - b6*(D4L_MB{-1}-ss_D4L_MB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%output persistence;
%b1 varies between 0.1 (extremely flexible) and 0.95(extremely persistent)
s.b1 = 0.3602;

%policy passthrough (the impact of monetary policy on real economy); 
%b2 varies between 0.1 (low impact) to 0.5 (strong impact)
s.b2 =0.1719;

%the impact of external demand on domestic output; 
%b3 varies between 0.1 and 0.7
s.b3 = 0.5996;

%the weight of the real interest rate and real exchange rate gaps in MCI;
%b4 varies from 0.3 to 0.8
s.b4 = 0.4556;

%Direct effect of remmittances on output gap (after controlling by foreign GDP) 
s.b5 = 0.15;

%Effect of money growth in MCI
s.b6 = 0.4536;


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
s.a1 = 0.4134;

% policy passthrough (the impact of rmc on inflation); 
% a2 varies between 0.1 (a flat Phillips curve and a high sacrifice ratio) 
% to 0.5 (a steep Phillips curve and a low sacrifice ratio)
s.a2 = 0.0434;

% the ratio of imported goods in firms' marginal costs (1-a3); 
% a3 varies between 0.9 for a closed economy to 0.5 for an open economy
s.a3 = 0.7739;

% Effect of imports/exports prices on non-core inflation
s.a4 = 0.1097; %%%%

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
s.e1 = 0.3894;

% Setting h2 equal to 0 implies float based on the UIP, h2 equal to 1 means
% managed FX following target appreciation.
s.h2 = 0.8167;

% D. Exchange Rate Target
% DLA_S_TAR = f1*DLA_S_TAR{-1} + (1-f1)*(D4L_CPI_TAR - ss_DLA_CPI_RW + ...
%             DLA_Z_BAR + f2*(D4L_CPI-D4L_CPI_TAR) + f3*L_GDP_GAP) + SHK_DLA_S_TAR;
% L_S_TAR = L_S_TAR{-1} + DLA_S_TAR/4;
s.f1 = 0.9723;
s.f2 = 0;
s.f3 = 0;

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
s.g1 = 0.9216;

% policy reactiveness: the weight put on inflation by the policy maker); 
% g2 has no upper limit but must be always higher then 0 (the Taylor principle)
s.g2 = 0.4556;

% policy reactiveness: the weight put on the output gap by the policy maker); 
% g3 has no upper limit but must be always higher then 0
s.g3 = 0.4503;

% degree to which the central bank does not control domestic money market
s.h1 = 0.0;

s.h3 = 0.2023;
% reaction to exchange rate devaluation with respect its target
s.g4 = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Money demand equation for monetary base or other monetary aggregate'
%D4L_MB = D4L_CPI + j1*D4L_GDP - j2*(RS-RS{-4}) - D4L_VEL;
s.j1 = 1.0;
s.j2 = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Speed of convergence of selected variables to their trend values.
% Used for risk premium, trends, foreign variables and commodity prices, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% persistent shock to risk premium
% SHKN_PREM = rho_SHKN_PREM*SHKN_PREM{-1} + SHK_PREM;
s.rho_SHKN_PREM         = 0.5;

% persistence in convergence of trend variables to their steady-state levels
% applies for:   DLA_GDP_BAR, DLA_Z_BAR, RR_BAR and RR_RW_BAR

s.rho_DLA_Z_BAR         = 0.9098;
s.rho_DLA_GDP_BAR       = 0.8178;
s.rho_RR_BAR            = 0.9103;
s.rho_RR_RW_BAR         = 0.8180;

% persistence in foreign GDP 
% L_GDP_RW_GAP = h2*L_GDP_RW_GAP{-1} + SHK_L_GDP_RW_GAP;
s.rho_L_GDP_RW_GAP      = 0.8037;
% adding persistence of remittance shock 
s.rho_REM_GDP           = 0.9884;
% adding persistence shock to D4L_CPI_NOSUBY 
s.rho_D4L_CPI_NOSUBY    = 0.6676; % 0.7998;
%adding persistence shock to money velocity
s.rho_PM_D4L_MB         = 0.7279;

% persistence in foreign interest rates (and inflation);
%RS_RW = rho_RS_RW*RS_RW{-1} + (1-rho_RS_RW)*(RR_BAR + DLA_CPI_RW) + SHK_RS_RW;
s.rho_RS_RW             = 0.8759;
s.rho_DLA_CPI_RW        = 0.6529;

% Speed of inflation target adjustment to the medium-term target (higher values mean slower adjustment)
% D4L_CPI_TAR = f1*D4L_CPI_TAR{-1} + (1-f1)*ss_D4L_CPI_TAR + SHK_D4L_CPI_TAR;
s.rho_D4L_CPI_TAR       = 0.8723;

% Imports/exports AR(2) dynamics
s.rho_D4L_IPEI_1        = 1.195;
s.rho_D4L_IPEI_2        = -0.5;

% Persistencia de variables de ajuste (baja persistencia para que converjan
% a 0 rápido)
% Para variables de Ajuste
s.rho_D4L_CPI_ADJ = 0.05;
s.rho_RS_ADJ = 0.05;
s.rho_D4S4L_GDP_ADJ = 0.05;
s.rho_D4L_S_ADJ = 0.05;
 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibrated Steady States
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Foreign trend inflation or inflation target
s.ss_DLA_CPI_RW = 2;
% Trend level of domestic real interest rate 
s.ss_RR_BAR = 1.0; %0.5; 
% Trend change in the real ER (negative number = real appreciation)
s.ss_DLA_Z_BAR = -2.0; 
% Potential output growth
s.ss_DLA_GDP_BAR = 3.5; %1
% Trend level of foreign real interest rate
s.ss_RR_RW_BAR = 0.50; %0.5;
% Domestic inflation target
s.ss_D4L_CPI_TAR = 4;  
% Steady state value for Remittance to GDP ratio
s.ss_REM_GDP  = 15.0;
%adding money growth steady state
s.ss_D4L_MB = 7.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard Deviation for Shocks in QPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s.std_SHK_L_GDP_GAP         = 1.165;
s.std_SHK_DLA_GDP_BAR       = 0.6435;
s.std_SHK_D4L_CPI_TAR       = 0; 
s.std_SHK_DLA_CPIXFE        = 1.0792;
s.std_SHK_D4L_CPI_NOSUBY    = 0.947;
s.std_SHK_D4L_MB            = 4.6364;
s.std_SHK_L_S               = 1.1758;
s.std_SHK_DLA_S_TAR         = 0.1073;
s.std_SHK_DLA_Z_BAR         = 0.8326;
s.std_SHK_PREM              = 1;
s.std_SHK_RS                = 0.4157;
s.std_SHK_RR_BAR            = 1.1431;
s.std_SHK_REM_GDP           = 0.536;
s.std_SHK_REM_GDP_BAR       = 0.2 * s.std_SHK_REM_GDP;
s.std_SHK_L_GDP_RW_GAP      = 1.3236;
s.std_SHK_RS_RW             = 0.2796;
s.std_SHK_DLA_CPI_RW        = 1.0357;
s.std_SHK_RR_RW_BAR         = 1.7438;
s.std_SHK_D4L_IPEI          = 3.36;
% Variables de ajuste
s.std_SHK_D4L_CPI_ADJ       = 0.01;
s.std_SHK_RS_ADJ            = 0.01;
s.std_SHK_D4S4L_GDP_ADJ     = 0.01;
s.std_SHK_D4L_S_ADJ         = 0.01;