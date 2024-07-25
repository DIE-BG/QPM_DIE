s = struct();
%{
Parametrización QPM - Versión 2 DIE

DIE
MJGM 07/24
%}

%{ 
% 1. Aggregate demand equation (the IS curve)
% L_GDP_GAP = b1*L_GDP_GAP{-1} - b2*MCI + b3*L_GDP_RW_GAP + ...
%     b5*(REM_GDP-ss_REM_GDP) + SHK_L_GDP_GAP;

% Real monetary conditions index (mci)
MCI = b4*(RR_GAP ) + b6*(-L_Z_GAP) - (1-b4-b6)*(PM_D4L_MB);

%}

s.b1 = 0.3602;
s.b2 =0.1719;
s.b3 = 0.5996;
s.b4 = 0.4556;
s.b5 = 0.15;
s.b6 = 0.4536;

%{ 
 2. Aggregate supply equation (the Phillips curve)
 DLA_CPIXFE =  a1*DLA_CPIXFE{-1} + (1-a1)*(DLA_CPIXFE{1}) + a2*RMC + ...
               SHK_DLA_CPIXFE;

 Real marginal cost (rmc)
 RMC = a3*L_GDP_GAP + (1-a3)*(L_ZIMP_GAP);
%}

s.a1 = 0.4134;
s.a2 = 0.0434;
s.a3 = 0.7739;

%{
% 3. Uncovered Interest Rate Parity (UIP)
L_S = h2*(L_S{-1} + DLA_S_TAR/4) + ...
      (1-h2)*(...
              (1-e1)*L_S{+1} +...
              e1*(L_S{-1} + 2/4*(D4L_CPI_TAR - ss_DLA_CPI_RW + DLA_Z_BAR)) + ...
              (RS_RW - RS + h3*(PM_D4L_MB) + PREM)/4)+...
     SHK_L_S;

DLA_S_TAR = f1*DLA_S_TAR{-1} + ...
            (1-f1)*(D4L_CPI_TAR - ss_DLA_CPI_RW + DLA_Z_BAR - ...
                    f2*(D4L_CPI-D4L_CPI_TAR) - f3*L_GDP_GAP)...
             + SHK_DLA_S_TAR;
%}

s.e1 = 0.3894;
s.h2 = 0.8167;
s.f1 = 0.9723;
s.f2 = 0;
s.f3 = 0;

%{ 
4. Monetary policy reaction function 
RS = h1*(4*(L_S{+1} - L_S) + RS_RW + PREM) + (1 - h1)*(g1*RS{-1} + ...
             (1 - g1)*(RSNEUTRAL + g2*(D4L_CPI{+4} - D4L_CPI_TAR{+4}) + ...
             g3*L_GDP_GAP + g4*(DLA_S{-1}-DLA_S_TAR{-1}))) + SHK_RS;
%}
% policy persistence; 
s.g1 = 0.9216;
% policy reactiveness: the weight put on inflation by the policy maker); 
% g2 has no upper limit but must be always higher then 0 (the Taylor principle)
s.g2 = 0.4556;
% policy reactiveness: the weight put on the output gap by the policy maker); 
% g3 has no upper limit but must be always higher then 0
s.g3 = 0.4503;

s.h1 = 0.0;
s.h3 = 0.2023;
s.g4 = 0; 

% Money demand equation for monetary base or other monetary aggregate'
% D4L_MB = D4L_CPI + j1*D4L_GDP - j2*(RS-RS{-4}) - D4L_VEL;
s.j1 = 1.0;
s.j2 = 0.1;

%% Procesos autorregresivos
s.rho_SHKN_PREM         = 0.5;
s.rho_DLA_Z_BAR         = 0.9098;
s.rho_DLA_GDP_BAR       = 0.8178;
s.rho_RR_BAR            = 0.9103;
s.rho_RR_RW_BAR         = 0.8180;
s.rho_L_GDP_RW_GAP      = 0.8037;
s.rho_REM_GDP           = 0.9884;
s.rho_D4L_CPI_NOSUBY    = 0.7998;
s.rho_PM_D4L_MB         = 0.7279;
s.rho_RS_RW             = 0.8759;
s.rho_DLA_CPI_RW        = 0.6529;
s.rho_D4L_CPI_TAR       = 0.8723;
% Para variables de Ajuste
s.rho_i_adj = 0.05;
s.rho_D4_GDP_SM_ADJ = 0.05;
s.rho_D4L_S_ADJ = 0.05;
%% Calibrated Steady States
s.ss_DLA_CPI_RW = 2;
s.ss_RR_BAR = 1.0; %0.5; 
s.ss_DLA_Z_BAR = -2.0; 
s.ss_DLA_GDP_BAR = 3.5; %1
s.ss_RR_RW_BAR = 0.50; %0.5;
s.ss_D4L_CPI_TAR = 4;  
s.ss_REM_GDP  = 15.0;
s.ss_D4L_MB = 7.5;

%% Standard Deviation for Shocks in QPM

s.std_SHK_L_GDP_GAP         = 1.165;
s.std_SHK_DLA_GDP_BAR       = 0.6435;
s.std_SHK_D4L_CPI_TAR       = 0; 
s.std_SHK_DLA_CPIXFE        = 1.0792;
s.std_SHK_D4L_CPI_NOSUBY    = 1.0452;
s.std_SHK_D4L_MB            = 4.6364;
s.std_SHK_L_S               = 1.1758;
s.std_SHK_DLA_S_TAR         = 0.1073;
s.std_SHK_DLA_Z_BAR         = 0.8326;
s.std_SHK_PREM              = 1;
s.std_SHK_RS                = 0.4157;
s.std_SHK_RR_BAR            = 1.1431;
s.std_SHK_REM_GDP           = 0.536;
s.std_SHK_L_GDP_RW_GAP      = 1.3236;
s.std_SHK_RS_RW             = 0.2796;
s.std_SHK_DLA_CPI_RW        = 1.0357;
s.std_SHK_RR_RW_BAR         = 1.7438;
