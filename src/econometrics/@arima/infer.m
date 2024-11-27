function [varargout] = infer(Mdl, Y, varargin)
%INFER Infer ARIMA model innovations and conditional variances
%
% Syntax:
%
%   [E,V,logL] = infer(Mdl,Y)
%   [E,V,logL] = infer(Mdl,Y,param1,val1,param2,val2)
%
% Description:
%
%   Infer the residuals and conditional variances of a univariate time 
%   series whose structure is characterized by an ARIMA model.
%
% Input Arguments:
%
%   Mdl - ARIMA model specification object, as produced by the ARIMA 
%     constructor or ARIMA/ESTIMATE method.
%
%   Y - Response data whose residuals and variances are inferred. Y is a
%     numObs-by-numPaths matrix, and represents the time series 
%     characterized by the model specification Mdl and is the continuation 
%     of the presample series Y0 (see below). Y is a column vector or a 
%     matrix. As a column vector, Y represents a single path of the 
%     underlying series; as a matrix, Y represents numObs observations of 
%     numPaths paths of an underlying time series, with observations across 
%     any row assumed to occur at the same time. The last observation of 
%     any series is assumed to be the most recent.
%
% Optional Input Parameter Name/Value Pairs:
%
%   'Y0'         Presample response data, providing initial values for the 
%                model. Y0 is a column vector or a matrix. If Y0 is a column 
%                vector, then it is applied to each inferred path. If Y0 is 
%                a matrix, then it must have at least numPaths columns. Y0 
%                may have any number of rows, provided at least Mdl.P 
%                observations exist to initialize the model. If the number 
%                of rows exceeds Mdl.P, then only the most recent Mdl.P
%                observations are used. If the number of columns exceeds 
%                numPaths, then only the first numPaths columns are used. 
%                If Y0 is unspecified, any necessary observations are 
%                backcasted (i.e., backward forecasted). The last row contains
%                the most recent observation.
%
%   'E0'         Mean-zero presample innovations, providing initial values 
%                for the model. E0 is a column vector or a matrix. If E0 is 
%                a column vector, then it is applied to each inferred path. 
%                If E0 is a matrix, then it must have at least numPaths 
%                columns. E0 may have any number of rows, provided sufficient
%                observations exist to initialize the ARIMA model as well 
%                as any conditional variance model (the number of observations
%                required is at least Mdl.Q, but may be more if a conditional 
%                variance model is included). If the number of rows exceeds 
%                the number necessary, then only the most recent observations
%                are used. If the number of columns exceeds numPaths, then 
%                only the first numPaths columns are used. If E0 is 
%                unspecified, any necessary observations are set to zero. 
%                The last row contains the most recent observation.
%
%   'V0'         Positive presample conditional variances, providing initial
%                values for any conditional variance model; if the variance 
%                of the model is constant, then V0 is unnecessary. V0 is a 
%                column vector or a matrix. If V0 is a column vector, then 
%                it is applied to each inferred path. If V0 is a matrix, 
%                then it must have at least numPaths columns. V0 may have 
%                any number of rows, provided sufficient observations exist 
%                to initialize the variance model. If the number of rows 
%                exceeds the number necessary, then only the most recent 
%                observations are used. If the number of columns exceeds 
%                numPaths, then only the first numPaths columns are used. 
%                If V0 is unspecified, any necessary observations are set 
%                to the unconditional variance of the conditional variance 
%                process. The last row contains the most recent observation.
%
%   'X'          Matrix of predictor data used to include a regression 
%                component in the conditional mean. Each column of X is a
%                separate time series, and the last row of each contains 
%                the most recent observation of each series. When presample 
%                responses Y0 are specified, the number of observations in 
%                X must equal or exceed the number of observations in Y; in 
%                the absence of presample responses, the number of observations
%                in X must equal or exceed the number of observations in Y 
%                plus Mdl.P. When the number of observations in X exceeds 
%                the number necessary, only the most recent observations 
%                are used. If missing, the conditional mean will have no 
%                regression component regardless of the presence of any 
%                regression coefficients found in the model.
%
% Output Arguments:
%
%   E - numObs-by-numPaths matrix of residuals inferred from the input 
%     series Y.
%
%   V - numObs-by-numPaths matrix of conditional variances inferred from
%     the input series Y.
%
%   logL - numPaths element vector of loglikelihood objective function 
%     values associated with the model specification Mdl. Each element of 
%     logL is associated with the corresponding path of Y.
%
% Notes:
%
%   o Missing values, indicated by NaNs, are removed from Y and X by listwise 
%     deletion (i.e., Y and X are merged into a composite series, and any row
%     of the combined series with at least one NaN is removed), reducing the 
%     effective sample size. Similarly, missing values in the presample data 
%     Y0, E0, and V0 are also removed by listwise deletion (Y0, E0, and V0 
%     are merged into a composite series, and any row of the combined series 
%     with at least one NaN is removed). The Y and X series, as well as the 
%     presample data, are also synchronized such that the last (most recent) 
%     observation of each component series occurs at the same time.
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
% References:
%
%   [1] Bollerslev, T. "Generalized Autoregressive Conditional 
%       Heteroskedasticity." Journal of Econometrics. Vol. 31, 1986, pp.
%       307-327.
%
%   [2] Bollerslev, T. "A Conditionally Heteroskedastic Time Series Model
%       for Speculative Prices and Rates of Return." The Review Economics
%       and Statistics. Vol. 69, 1987, pp 542-547.
%
%   [3] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
%       Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
%       NJ: Prentice-Hall, 1994.
%
%   [4] Enders, W. Applied Econometric Time Series. Hoboken, NJ: John Wiley
%       & Sons, 1995.
%
%   [5] Engle, R. F. "Autoregressive Conditional Heteroskedasticity with
%       Estimates of the Variance of United Kingdom Inflation." 
%       Econometrica. Vol. 50, 1982, pp. 987-1007.
%
%   [6] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
% See also ARIMA, ESTIMATE, FORECAST, SIMULATE.

% Copyright 2018 The MathWorks, Inc.   

if numel(Mdl) > 1
   error(message('econ:arima:infer:NonScalarModel'))
end

%
% Check input parameters and set defaults.
%

if (nargin < 2) || isempty(Y)
   error(message('econ:arima:infer:InvalidSeries'))
end

parser = inputParser;
parser.addRequired ('requiredY',    @(x) validateattributes(x, {'double'}, {}, '', 'response data'));
parser.addParameter('Y0'       , 0, @(x) validateattributes(x, {'double'}, {}, '', 'presample responses'));
parser.addParameter('E0'       , 0, @(x) validateattributes(x, {'double'}, {}, '', 'presample residuals'));
parser.addParameter('V0'       , 0, @(x) validateattributes(x, {'double'}, {}, '', 'presample variances'));
parser.addParameter('X'        ,[], @(x) validateattributes(x, {'double'}, {}, '', 'regression matrix'));

try
   parser.parse(Y, varargin{:});
catch exception
  exception.throwAsCaller();
end

Y  = parser.Results.requiredY;
Y0 = parser.Results.Y0;
E0 = parser.Results.E0;
V0 = parser.Results.V0;
X  = parser.Results.X;

%
% Get model parameters and extract lags associated with non-zero coefficients. 
%
% In the code segment below, AR and MA represent the compound auto-regressive 
% and moving average polynomials, respectively, including the effects of 
% integration and seasonality.
%

constant              = Mdl.PrivateConstant;    % Additive constant
variance              = Mdl.PrivateVariance;    % Variance
isVarianceConditional = any(strcmp(class(variance), {'garch' 'gjr' 'egarch'}));  % Is it a conditional variance model (e.g., GARCH)?

if isVarianceConditional                        % Allow for a conditional variance model
   P = max(Mdl.PrivateVariance.P, Mdl.P);
   Q = max(Mdl.PrivateVariance.Q, Mdl.Q);       % Total number of lagged e(t) needed
else
   P = Mdl.P;                                   % Total number of lagged y(t) needed
   Q = Mdl.Q;                                   % Total number of lagged e(t) needed
end

AR     = getLagOp(Mdl, 'Compound AR');
MA     = reflect(getLagOp(Mdl, 'Compound MA')); % This negates the MA coefficients
LagsAR = AR.Lags;                               % Lags of non-zero AR coefficients
LagsMA = MA.Lags;                               % Lags of non-zero MA coefficients
LagsMA = LagsMA(LagsMA > 0);                    % Exclude lag zero (MA terms only)

if isempty(LagsAR)
   AR = [];
else
   AR = AR.Coefficients;                        % Lag Indexed Array
   AR = [AR{LagsAR}];                           % Non-zero AR coefficients (vector)
end

if isempty(LagsMA)
   MA = [];
else
   MA  = MA.Coefficients;                       % Lag Indexed Array
   MA  = [MA{LagsMA}];                          % Non-zero MA coefficients (vector)
end

%
% Ensure coefficients are specified.
%

if any(isnan(constant))
   error(message('econ:arima:infer:UnspecifiedConstant'))
end

if any(isnan(AR))
   error(message('econ:arima:infer:UnspecifiedAR'))
end

if any(isnan(MA))
   error(message('econ:arima:infer:UnspecifiedMA'))
end

isRegressionIncluded = ~any(strcmpi('X', parser.UsingDefaults));

if isRegressionIncluded
   beta = Mdl.PrivateBeta;
   if isempty(beta) || any(isnan(beta))
      error(message('econ:arima:infer:UnspecifiedBeta'))
   end
   if numel(beta) ~= size(X,2)
      error(message('econ:arima:infer:InconsistentRegression'))
   end
end

if ~isVarianceConditional && any(isnan(variance))
   error(message('econ:arima:infer:UnspecifiedVariance'))
end

%
% Remove missing observations (NaN's) via listwise deletion.
%

if any(isnan(Y0(:))) || any(isnan(E0(:))) || any(isnan(V0(:)))
   [Y0, E0, V0] = internal.econ.LagIndexableTimeSeries.listwiseDelete(Y0, E0, V0);  % Pre-sample data
end

if any(isnan(Y(:))) || any(isnan(X(:)))
   [Y,X] = internal.econ.LagIndexableTimeSeries.listwiseDelete(Y, X);               % In-sample data
end

%
% Compute the total number of observations generated for each path as the
% sum of the number of observations of the input series y(t) and the number 
% of presample observations needed to initialize the recursions.
%

[numObs, numPaths] = size(Y);
maxPQ    = max([P Q]);              % Maximum presample lags needed   
T        = numObs + maxPQ;          % Total number of observations generated

%
% Check any user-specified regression data for sufficient observations.
%

isY0specified = ~any(strcmpi('Y0', parser.UsingDefaults));

if isRegressionIncluded

   try
     if isY0specified
%
%       Presample responses (Y0) are specified, and so the number of observations 
%       in X must equal or exceed the number of observations in Y.
%
        X = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(T,size(X,2)), 'X', X, numObs);
     else
%
%       No presample responses (Y0) are specified, and so the number of 
%       observations in X must equal or exceed the number of observations 
%       in Y plus the degree of the autoregressive polynomial (Mdl.P).
%
        X = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(T,size(X,2)), 'X', X, numObs + Mdl.P);
     end

   catch exception
     error(message('econ:arima:infer:InsufficientXRows', numObs + (Mdl.P * (~isY0specified))))
   end

end

%
% Check any user-specified presample observations used for conditioning, or 
% generate any required observations automatically.
%

isE0specified = ~any(strcmpi('E0', parser.UsingDefaults));

if isE0specified      % Did the user specify presample e(t) observations?

%
%  Check user-specified presample data for the residuals e(t).
%

   E0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(maxPQ,numPaths), 'E0', E0, Q);

else

%
%  The user did not specify presample e(t) observations, so initialize any 
%  required presample observations with the unconditional mean of zero.
%

   E0 = zeros(maxPQ,numPaths);     % Unconditional mean of e(t)

end

if isY0specified                   % Did the user specify presample y(t) observations?

%
%  Check user-specified presample data for the responses y(t).
%

   Y0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(ones(maxPQ,numPaths), 'Y0', Y0, Mdl.P);

else

%
%  The user did not specify presample y(t) observations, so backcast the 
%  initial values. 
%

   if maxPQ > 0

      if isRegressionIncluded      % ARIMAX model
%
%        Since Y(t) and X(t) are synchronized such that the most recent
%        observations coincide, the first max(P,Q) observations of X(t) - when 
%        time-reversed - serve as the forecast of X(t) (i.e., XF) into the past.
%
         Y0 = forecast(Mdl, maxPQ, 'Y0', Y(end:-1:1,:)          , 'E0', zeros(maxPQ,numPaths), ...
                                   'X0', X(end:-1:(1 + maxPQ),:), 'XF', X(maxPQ:-1:1,:));

      else                         % ARIMA model (no regression component)

         Y0 = forecast(Mdl, maxPQ, 'Y0', Y(end:-1:1,:), 'E0', zeros(maxPQ,numPaths)); 
      end

      Y0 = Y0(end:-1:1,:);

   else

      Y0 = [];
   end

end

%
% Size the responses y(t) and residuals e(t), initialize each with the 
% presample data, and transpose for efficiency.
%

Y = [Y0'  Y'];
E = [E0'  zeros(numPaths,numObs)];
X = X';

%
% Infer the residuals and conditional variances, and calculate values of
% the log-likelihood function.
%

coefficients = [-constant  AR  MA]';

if isVarianceConditional                % Conditional variance model

%
%  Infer the residuals.
%

   if isRegressionIncluded              % ARIMAX model

      E = internal.econ.arimaxMex(LagsAR, LagsMA, coefficients, maxPQ, Y', E', X', -beta)';
%
%     Note:
%
%     To improve runtime performance, the following code segment has been replaced 
%     by the private/undocumented C-MEX file called above, which does exactly the 
%     same thing. To restore the original MATLAB code, simply uncomment the 
%     following code segment:

%       I = ones(numPaths,1);
%
%       for t = (maxPQ + 1):T
%           data   = [I  Y(:,t - LagsAR)  E(:,t - LagsMA)];
%           E(:,t) = data * coefficients  -  beta * X(:,t);
%       end
      
   else                                 % ARIMA model (no regression component)

      E = internal.econ.arimaxMex(LagsAR, LagsMA, coefficients, maxPQ, Y', E')';
%
%     Note:
%
%     To improve runtime performance, the following code segment has been replaced 
%     by the private/undocumented C-MEX file called above, which does exactly the 
%     same thing. To restore the original MATLAB code, simply uncomment the 
%     following code segment:

%       I = ones(numPaths,1);
%
%       for t = (maxPQ + 1):T
%           data   = [I  Y(:,t - LagsAR)  E(:,t - LagsMA)];
%           E(:,t) = data * coefficients;
%       end

   end
   
   varargout(1) = {E(:,(maxPQ + 1):T)'};

   if nargout > 1
%
%     Infer the conditional variances and log-likelihood values.
%
      isV0specified = ~any(strcmpi('V0', parser.UsingDefaults));

      if isE0specified && isV0specified
         [varargout{2:nargout}] = infer(variance, E(:,(maxPQ + 1):T)', 'E0', E0, 'V0', V0);
      elseif isE0specified
         [varargout{2:nargout}] = infer(variance, E(:,(maxPQ + 1):T)', 'E0', E0);
      elseif isV0specified
         [varargout{2:nargout}] = infer(variance, E(:,(maxPQ + 1):T)', 'V0', V0);
      else
         [varargout{2:nargout}] = infer(variance, E(:,(maxPQ + 1):T)');
      end

   end

else                              % Constant variance model

   V = variance(ones(numPaths,T));
   
   if isRegressionIncluded
      coefficients = [coefficients ; beta'];
   end

   isDistributionT = strcmpi(Mdl.PrivateDistribution.Name, 'T'); % Is the distribution Student's t?

   if isDistributionT && isnan(Mdl.PrivateDistribution.DoF) && (nargout > 2)
%
%     Throw an error only if necessary. Even for t distributions, the degrees 
%     of freedom parameter is only needed to evaluate the log-likelihood, and 
%     not to infer the residuals.
%
      error(message('econ:arima:infer:UnspecifiedDoF'))
   end

   if isDistributionT
      [nLogL,~,E] = arima.nLogLikeT([coefficients ; Mdl.PrivateDistribution.DoF], Y, X, E, V, LagsAR, LagsMA, numPaths, maxPQ, T);
   else
      [nLogL,~,E] = arima.nLogLikeGaussian(coefficients, Y, X, E, V, LagsAR, LagsMA, numPaths, maxPQ, T);
   end

%
%  Strip the initialization period for output purposes.
%

   varargout = {E(:,(maxPQ + 1):T)'  V(:,(maxPQ + 1):T)'  -nLogL'};
end

end