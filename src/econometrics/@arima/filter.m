function [Y,E,V] = filter(Mdl, Z, varargin)
%FILTER Filter disturbances through an ARIMA model
%
% Syntax:
%
%   [Y,E,V] = filter(Mdl,Z)
%   [Y,E,V] = filter(Mdl,Z,param1,val1,param2,val2,...)
%
% Description:
%
%   Filter user-specified disturbances to produce responses, innovations, 
%   and conditional variances of a univariate ARIMA model.
%
% Input Arguments:
%
%   Mdl - ARIMA model specification object, as produced by the ARIMA 
%     constructor or ARIMA/ESTIMATE method.
%
%   Z - A numObs-by-numPaths matrix of disturbances used to drive the output
%     innovations process E(t) such that, for a variance process V(t), 
%     E(t) = sqrt(V(t))*Z(t). Z represents the continuation of the presample
%     series Z0 (see below). Z is a column vector or a matrix. As a column 
%     vector, Z represents a single path of the underlying disturbance series; 
%     as a matrix, Z represents numObs observations of numPaths paths of the 
%     underlying disturbance series, with observations across any row assumed 
%     to occur at the same time. The last observation of any series is 
%     assumed to be the most recent.
%
% Optional Input Parameter Name/Value Pairs:
%
%   'Y0'         Presample response data, providing initial values for the 
%                model. Y0 is a column vector or a matrix. As a column 
%                vector, Y0 is applied to each output path. As a matrix, Y0
%                must have at least numPaths columns. Y0 may have any number 
%                of rows, provided at least Mdl.P observations exist to 
%                initialize the model. If the number of rows exceeds Mdl.P, 
%                then only the most recent Mdl.P observations are used. If 
%                the number of columns exceeds numPaths, then only the first 
%                numPaths columns are used. If Y0 is unspecified, any 
%                necessary presample observations are set to the unconditional
%                mean for stationary AR processes, and to zero if the process 
%                is non-stationary or contains a regression component. The 
%                last row contains the most recent observation.
%
%   'Z0'         Presample disturbances providing initial values for the 
%                input disturbance series Z (see above). Z0 is a column 
%                vector or a matrix. As a column vector, Z0 is applied to
%                each output path. As a matrix, Z0 must have at least 
%                numPaths columns. Z0 may have any number of rows, provided 
%                sufficient observations exist to initialize the moving 
%                average component of the ARIMA model as well as any 
%                conditional variance model (the number of observations 
%                required is at least Mdl.Q, but may be more if a conditional 
%                variance model is included). If the number of rows exceeds 
%                the number necessary, then only the most recent observations 
%                are used. If the number of columns exceeds numPaths, then 
%                only the first numPaths columns are used. If no presample 
%                data is specified, any necessary observations are set to 
%                zero. The last row contains the most recent observation.
%
%   'V0'         Positive presample conditional variances, providing initial
%                values for the model. V0 is a column vector or a matrix. As
%                a column vector, V0 is applied to each simulated path. As 
%                a matrix, V0 must have at least numPaths columns. V0 may 
%                have any number of rows, provided sufficient observations 
%                exist to initialize the moving average component of the 
%                ARIMA model as well as any conditional variance model (the
%                number of observations required is at least Mdl.Q, but
%                may be more if a conditional variance model is included).
%                If the number of rows exceeds the number necessary, then 
%                only the most recent observations are used. If the number 
%                of columns exceeds numPaths, then only the first numPaths 
%                columns are used. If no presample variance data is specified,
%                any necessary observations are set to the unconditional 
%                variance of the conditional variance process. The last row 
%                contains the most recent observation.
%
%   'X'          Matrix of predictor data used to include a regression 
%                component in the conditional mean. Each column of X is a
%                separate time series, and the last row of each contains
%                the most recent observation of each series. The number of
%                observations in X must equal or exceed the number of 
%                observations in Z. When the number of observations in X 
%                exceeds the number necessary, only the most recent 
%                observations are used. If missing, the conditional mean 
%                will have no regression component regardless of the presence 
%                of any regression coefficients found in the model.
%
% Output Arguments:
%
%   Y - numObs-by-numPaths matrix of simulated responses, and the continuation
%     of the presample series Y0 (see above).
%
%   E - numObs-by-numPaths matrix of simulated innovations with conditional 
%     variances V (see below). E is a scaled innovation, or disturbance, 
%     series such that E(t) = sqrt(V(t))*Z(t).
%
%   V - numObs-by-numPaths matrix of conditional variances of the innovations 
%     in E such that E(t) = sqrt(V(t))*Z(t), and the continuation of the 
%     presample series V0 (see above).
%
% Notes:
%
%   o The FILTER method is designed to generalize the SIMULATE method, and 
%     both filter a series of disturbances to produce output responses, 
%     innovations, and conditional variances. However, whereas SIMULATE 
%     auto-generates a series of mean-zero, unit-variance, independent and 
%     identically-distributed (iid) disturbances Z(t) according to the 
%     distribution found in the model Mdl, FILTER allows users to directly 
%     specify their own Z(t) disturbances.
%
%   o Missing values, indicated by NaNs, are removed from Z and X by listwise 
%     deletion (i.e., Z and X are merged into a composite series, and any row
%     of the combined series with at least one NaN is removed), reducing the 
%     effective sample size. Missing values in the presample data Y0, Z0, and
%     V0 are also removed by listwise deletion. Y0, Z0, and V0 are merged 
%     into a composite series, and any row of the combined series with at 
%     least one NaN is removed. The presample data is also synchronized such 
%     that the last (most recent) observation of each series occurs at the 
%     same time.
%
%  o  Regression models included in the conditional mean are based on the
%     presence of the predictor matrix X. Although each column of the output 
%     time series represents a different path of the corresponding univariate 
%     stochastic process, the regression matrix X represents as a single 
%     path of a (possibly) multivariate time series matrix in which each 
%     column is a different time series. When the conditional mean has a 
%     regression component, the entire predictor matrix X is applied to 
%     every column of the output time series. 
%
% Examples:
%
%   o Simulate-to-Filter Round-Trip:
%     Consider a "round-trip" example to illustrate the relationship between
%     SIMULATE and FILTER. An ARIMA(2,0,1) model is first simulated, then
%     the simulated model innovations (e) are first standardized and then 
%     filtered:
%
%     model   = arima('Constant',0,'AR',{0.5 -0.8},'MA',-0.5,'Variance',0.1);
%     [y,e,v] = simulate(model, 100);       % Simulate via Monte Carlo
%     Z       = e./sqrt(v);                 % Standardize the innovations
%     [Y,E,V] = filter(model, Z);           % Filter Z(t)
%
%     In the above round-trip example, the corresponding outputs of
%     SIMULATE and FILTER are identical.
%
%   o Impulse Response:
%     The impulse response assesses the dynamic behavior of a system to a 
%     one-time, unit impulse. The following example illustrates the first
%     20 responses of a mean-zero ARIMA(2,0,1) model to a unit shock:
%
%     model = arima('Constant',0,'AR',{0.5 -0.8},'MA',-0.5,'Variance',0.1);
%     Y     = filter(model, [1 ; zeros(19,1)], 'Y0', zeros(model.P,1));
%     Y     = Y / Y(1);     % Normalize to ensure that Y(1) = 1
%     stem((0:numel(Y)-1)', Y, 'filled'); title('Impulse Response')
%
%   o Step Response:
%     The step response assesses the dynamic behavior of a system to a 
%     persistent change in a variable. The following example illustrates 
%     the first 20 responses of a mean-zero ARIMA(2,0,1) model to a sequence
%     of unit disturbances:
%
%     model = arima('Constant',0,'AR',{0.5 -0.8},'MA',-0.5,'Variance',0.1);
%     Y     = filter(model, ones(20,1), 'Y0', zeros(model.P,1));
%     Y     = Y / Y(1);     % Normalize to ensure that Y(1) = 1
%     stem((0:numel(Y)-1)', Y, 'filled'); title('Step Response')
%
% References:
%
%   [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
%       Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
%       NJ: Prentice-Hall, 1994.
%
%   [2] Enders, W. Applied Econometric Time Series. Hoboken, NJ: John Wiley
%       & Sons, 1995.
%
%   [3] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
% See also ARIMA, SIMULATE.

% Copyright 2018 The MathWorks, Inc.

if numel(Mdl) > 1
   error(message('econ:arima:filter:NonScalarModel'))
end

%
% Check input parameters and set defaults.
%

if nargin < 2
   error(message('econ:arima:filter:NotEnoughInputs'))
end

parser = inputParser;
parser.addRequired ('requiredZ',    @(x) validateattributes(x, {'double'}, {'nonempty'}, '', 'standardized disturbances'));
parser.addParameter('Y0'       , 0, @(x) validateattributes(x, {'double'}, {}          , '', 'presample responses'));
parser.addParameter('Z0'       , 0, @(x) validateattributes(x, {'double'}, {}          , '', 'presample disturbances'));
parser.addParameter('V0'       , 0, @(x) validateattributes(x, {'double'}, {}          , '', 'presample variances'));
parser.addParameter('X'        ,[], @(x) validateattributes(x, {'double'}, {}          , '', 'regression matrix'));

try
  parser.parse(Z, varargin{:});
catch exception
  exception.throwAsCaller();
end

Z  = parser.Results.requiredZ;
Y0 = parser.Results.Y0;
Z0 = parser.Results.Z0;
V0 = parser.Results.V0;
X  = parser.Results.X;

%
% Get model parameters and extract lags associated with non-zero coefficients.
%
% In the code segment below, AR and MA represent the compound auto-regressive 
% and moving average polynomials, respectively, including the effects of 
% integration and seasonality.
%

constant = Mdl.PrivateConstant;       % Additive constant
variance = Mdl.PrivateVariance;       % Conditional variance

isVarianceConditional = any(strcmp(class(variance), {'garch' 'gjr' 'egarch'})); % Is it a conditional variance model?

if isVarianceConditional                   % Allow for a conditional variance model
   P = max(Mdl.PrivateVariance.P, Mdl.P);
   Q = max(Mdl.PrivateVariance.Q, Mdl.Q);  % Total number of lagged z(t) needed
else
   P = Mdl.P;                              % Total number of lagged y(t) needed
   Q = Mdl.Q;                              % Total number of lagged z(t) needed
end

AR     = Mdl.getLagOp('Compound AR');
MA     = Mdl.getLagOp('Compound MA');
LagsAR = AR.Lags;                     % Lags of non-zero AR coefficients
LagsAR = LagsAR(LagsAR > 0);
LagsMA = MA.Lags;                     % Lags of non-zero MA coefficients

isARstable = isStable(AR);            % Determine if the process is AR stable

if isempty(LagsAR)
   AR =  [];
else
   AR = AR.Coefficients;              % Lag Indexed Array
   AR = -[AR{LagsAR}];                % Non-zero AR coefficients (vector)
end

MA = MA.Coefficients;                 % Lag Indexed Array
MA = [MA{LagsMA}];                    % Non-zero MA coefficients (vector)

%
% Ensure coefficients are specified.
%

if any(isnan(constant))
   error(message('econ:arima:filter:UnspecifiedConstant'))
end

if any(isnan(AR))
   error(message('econ:arima:filter:UnspecifiedAR'))
end

if any(isnan(MA))
   error(message('econ:arima:filter:UnspecifiedMA'))
end

isRegressionIncluded = ~any(strcmpi('X', parser.UsingDefaults));

if isRegressionIncluded
   beta = Mdl.PrivateBeta;
   if isempty(beta) || any(isnan(beta))
      error(message('econ:arima:filter:UnspecifiedBeta'))
   end
   if numel(beta) ~= size(X,2)
      error(message('econ:arima:filter:InconsistentRegression'))
   end
end

if isVarianceConditional
   if isnan(variance.UnconditionalVariance)
      error(message('econ:arima:filter:UnspecifiedConditionalVariance'))
   end
else
   if any(isnan(variance))
      error(message('econ:arima:filter:UnspecifiedVariance'))
   end
end

%
% Remove missing observations (NaN's) via listwise deletion.
%

if any(isnan(Y0(:))) || any(isnan(Z0(:))) || any(isnan(V0(:)))
   [Y0, Z0, V0] = internal.econ.LagIndexableTimeSeries.listwiseDelete(Y0, Z0, V0); % Pre-sample data
end

if any(isnan(Z(:))) || any(isnan(X(:)))
   [Z,X] = internal.econ.LagIndexableTimeSeries.listwiseDelete(Z, X);              % In-sample data
end

%
% Compute the total number of observations generated for each path as the
% sum of the number of observations simulated and the number of presample 
% observations needed to initialize the recursion.
%

[numObs, numPaths] = size(Z);
maxPQ              = max([P Q]);     % Maximum presample lags needed   
T                  = numObs + maxPQ; % Total number of observations generated

%
% Check any user-specified presample observations used for conditioning, or 
% generate any required observations automatically.
%

isZ0specified = ~any(strcmpi('Z0', parser.UsingDefaults));

if isZ0specified      % Did the user specify presample z(t) observations?

%
%  Check user-specified presample data for the innovations z(t).
%

   Z0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(maxPQ,numPaths), 'Z0', Z0, Q);

else

%
%  The following check is included to avoid a forward compatibility issue
%  in the event presample z(t) observations are generated from the distribution 
%  found in the model Mdl (i.e., using the utility "simulateStandardizedVariates"
%  as for GARCH, EGARCH, and GJR classes).
%

   isDistributionT = strcmpi(Mdl.PrivateDistribution.Name, 'T'); % Is the distribution Student's t?

   if isDistributionT && isnan(Mdl.PrivateDistribution.DoF)
      error(message('econ:arima:filter:UnspecifiedDoF'))
   end

%
%  The user did not specify presample z(t) observations, so assume zero.
%

   Z0 = zeros(maxPQ,numPaths);

end

if ~any(strcmpi('Y0', parser.UsingDefaults))  % Did the user specify presample y(t) observations?

%
%  Check user-specified presample data for the responses y(t).
%

   Y0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(ones(maxPQ,numPaths), 'Y0', Y0, Mdl.P);

%
%  Size the responses y(t) and initialize with specified data.
%
   Y            = zeros(numPaths,T);
   Y(:,1:maxPQ) = Y0';

else

%
%  The user did not specify presample y(t) observations. 
%

   if isARstable && (sum(AR) ~= 1) && ~isRegressionIncluded
%
%     The model is AR-stable and without a regression component, so compute 
%     the unconditional (i.e., long-run) mean of the y(t) process directly 
%     from the parameters of the model and use it to initialize any required 
%     presample observations.
%
      average = constant / (1 - sum(AR));
      Y       = repmat([average(ones(1,maxPQ)) zeros(1,numObs)], numPaths, 1);

   else
%
%     The model is not AR-stable, and so a long-run mean of the y(t) process 
%     cannot be calculated from the model. The following simply assumes zeros 
%     for any required presample observations for y(t).
%
      Y  = zeros(numPaths,T);

   end

end

%
% Check any user-specified presample variances v(t), or auto-generate them.
%

isV0specified = ~any(strcmpi('V0', parser.UsingDefaults));

if isV0specified      % Did the user specify presample v(t) observations?

%
%  The following check is needed because the FILTER method requires initial
%  conditional variance observations (with or without a contained conditional 
%  variance model). That is, as opposed to SIMULATE, FILTER requires z(t), 
%  rather than e(t) = sqrt(v(t)) * z(t), and so FILTER requires V0 observations 
%  even in the absence of a conditional variance model.
%

   if isVarianceConditional   % Is there a conditional variance model?
%
%     In the presence of a contained conditional variance model, FILTER 
%     requires max(Mdl.Q, Mdl.Variance.P, Mdl.Variance.Q) initial variance 
%     observations.
%
      V0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(maxPQ,numPaths), 'V0', V0, max(Mdl.PrivateVariance.P,Q));
      
   else
%
%     In the absence of a contained conditional variance model, FILTER still
%     requires Mdl.Q variance observations to initialize the model.
%
      V0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(maxPQ,numPaths), 'V0', V0, Q);

   end

else

%
%  In the absence of user-specified presample variances, set any necessary
%  presample observations to the unconditional variance of any contained
%  variance model (consistent with the default approach of the contained
%  conditional variance model), or simply equal to the constant variance of
%  constant-variance models.
%
   if isVarianceConditional   % Is there a conditional variance model?
      V0 = repmat(variance.UnconditionalVariance, maxPQ, numPaths);
   else
      V0 = repmat(variance, maxPQ, numPaths);
   end

end

%
% Check any user-specified regression data for sufficient observations.
%

if isRegressionIncluded
   try
     X = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(T,size(X,2)), 'X', X, numObs);
     X = X.';
   catch exception
     error(message('econ:arima:filter:InsufficientXRows'))
   end
end

%
% Generate the innovations e(t).
%

E            = zeros(numPaths,T);
E(:,1:maxPQ) = (sqrt(V0) .* Z0)';

if isVarianceConditional   % Is there a conditional variance model?

   if isZ0specified && isV0specified 
      [V,e] = filter(variance, Z, 'Z0', Z0, 'V0', V0);
   elseif isZ0specified
      [V,e] = filter(variance, Z, 'Z0', Z0);
   elseif isV0specified
      [V,e] = filter(variance, Z, 'V0', V0);
   else
      [V,e] = filter(variance, Z);
   end

   E(:,(maxPQ + 1):end) = e';

else                       % Then it's a constant-variance model

%
%  Format the model innovations e(t) (i.e., the non-standardized disturbances)
%  such that e(t) = sqrt(v(t)) * z(t).
%

   E(:,(maxPQ + 1):end) = Z' * sqrt(variance);

   if nargout > 2
      V = variance(ones(numObs, numPaths));
   end

end

%
% Simulate the data y(t).
%

coefficients = [constant  AR  MA]';  % ARIMA coefficient vector
I            = ones(numPaths,1);

if isRegressionIncluded              % ARIMAX model

   for t = (maxPQ + 1):T
       data   = [I  Y(:,t - LagsAR)  E(:,t - LagsMA)];
       Y(:,t) = data * coefficients  +  beta * X(:,t);
   end

else                                 % ARIMA model (no regression component)

   for t = (maxPQ + 1):T
       data   = [I  Y(:,t - LagsAR)  E(:,t - LagsMA)];
       Y(:,t) = data * coefficients;
   end

end

%
% Since max(P,Q) observations have been prepended to the output processes 
% to compensate for presample effects, strip the start-up, retain only the 
% last numObs observations, and transpose to a conventional time series
% format.
%

Y = Y(:,(maxPQ + 1):T)';
E = E(:,(maxPQ + 1):T)';

end