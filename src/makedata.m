%% Load quarterly data
% Base de datos inicial
d = databank.fromCSV('data/data.csv');


%% Make log of variables
exceptions = {'RS','RS_RW'};

list = dbnames(d);

for i = 1:length(list)
    if isempty(strmatch(list{i},exceptions,'exact'))
        d.(['L_' list{i}]) = 100*log(d.(list{i}));
    end
end


%% Define the real exchange rate (en niveles)
d.L_Z = d.L_S + d.L_CPI_RW - d.L_CPIXFE;

%% Growth rate annualized qoq and yoy
exceptions = {'RS','RS_RW'};

list = dbnames(d);

for i = 1:length(list)
    if isempty(strmatch(list{i}, exceptions,'exact'))
        if length(list{i})>1
            % Condicional para Logs
            if strcmp('L_', list{i}(1:2))
                % Q-o-Q annualized growth rates
                d.(['DLA_' list{i}(3:end)])  = 4*(d.(list{i}) - d.(list{i}){-1});
                % Y-o-Y growth rates
                d.(['D4L_' list{i}(3:end)]) = d.(list{i}) - d.(list{i}){-4};
            end
        end
    end
end

%% Real Interest Rates
% Domestic real interest rate (history)
d.RR = d.RS - d.D4L_CPI;
% Foreign real interest rate (history)
d.RR_RW = d.RS_RW - d.D4L_CPI_RW;

%% Trends and Gaps - Hodrick-Prescott filter 
% For GDP_RW the observable for the QPM is L_GDP_RW_GAP

list = {'RR_RW','L_GDP_RW'};

for i = 1:length(list)
    [d.([list{i} '_BAR']), d.([list{i} '_GAP'])] = hpf(d.(list{i}));
end

%% Variables que podrían utilizarse 
% Trend and Gap for Output - Band-pass filter
% save HP results
% d.L_GDP_BAR_HP = d.L_GDP_BAR;
% d.L_GDP_GAP_HP = d.L_GDP_GAP; 
% 
% % Band-pass
% d.L_GDP_GAP = bpass(d.L_GDP,inf,[6,32],'detrend',false);
% d.L_GDP_BAR = hpf((d.L_GDP-d.L_GDP_GAP),inf,'lambda',5);
% d.DLA_GDP_BAR = 4*(d.L_GDP_BAR - d.L_GDP_BAR{-1});

% Growth rates of equilibria ('L_RWOIL_BAR', 'L_RWFOOD_BAR' excluded coherent with Model_MPAF4_est.mod)
% list = {'L_GDP_BAR', 'L_GDP_RW_BAR', 'L_Z_BAR', 'L_RWOIL_BAR', 'L_RWFOOD_BAR'};
% list = {'L_GDP_BAR', 'L_GDP_RW_BAR', 'L_Z_BAR'};
% 
% for i = 1:length(list)
%    d.(['DLA' list{i}(2:end)]) = 4*(d.(list{i}) - d.(list{i}){-1});
%    d.(['D4L' list{i}(2:end)]) = d.(list{i}) - d.(list{i}){-4};
% end

% Implied risk premium
% d.PREM = d.RR_BAR - d.RR_RW_BAR - d.DLA_Z_BAR;
% d.SHKN_PREM = tseries(get(d.L_S,'range'),0);


% Compute the exchange rate target over the history 
% Adjust lambda
% d.L_S_TAR     = hpf(d.L_S, Inf, 'lambda=', 1600);
% d.DLA_S_TAR   = 4*(d.L_S_TAR - d.L_S_TAR{-1});


%% Save the database
% Database is saved in file 'history.csv'
databank.toCSV(d,'output\history.csv', Inf);

