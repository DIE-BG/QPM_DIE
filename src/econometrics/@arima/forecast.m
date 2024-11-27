function [Y,YMSE,V] = forecast(Mdl, numPeriods, varargin)
%FORECAST Forecast ARIMA model responses and conditional variances
%
% Syntax:
%
%   [Y,YMSE,V] = forecast(Mdl,numPeriods,Y0)
%   [Y,YMSE,V] = forecast(Mdl,numPeriods,Y0,param1,val1,param2,val2,...)
%
% Description:
%
%   Forecast responses and conditional variances of a univariate time series
%   whose structure is characterized by an ARIMA model. 
%
% Input Arguments:
%
%   Mdl - ARIMA model specification object, as produced by the ARIMA 
%     constructor or ARIMA/ESTIMATE method.
%
%   numPeriods - Positive integer specifying the number of periods in the 
%     forecast horizon. 
%
%   Y0 - Presample response data, providing initial values for the model to 
%     forecast. Y0 is a column vector or a matrix. If Y0 is a column vector, 
%     then it is applied to each forecasted path. If Y0 is a matrix, then it 
%     must have numPaths columns (see notes below). Y0 may have any number 
%     of rows, provided at least Mdl.P observations exist to initialize the 
%     model. If the number of rows exceeds Mdl.P, then only the most recent 
%     Mdl.P observations are used. The last row contains the most recent 
%     observation.
%
% Optional Input Parameter Name/Value Pairs:
%
%   'E0'   Mean-zero presample innovations, providing initial values for the 
%          moving average (MA) component or conditional variance model. E0 
%          is a column vector or a matrix. If E0 is a column vector, then it 
%          is applied to each forecasted path. If E0 is a matrix, then it must
%          have numPaths columns (see notes below). E0 may have any number 
%          of rows, provided sufficient observations exist to initialize the 
%          MA component or any conditional variance model (the number of 
%          observations required is at least Mdl.Q, but may be more if a 
%          conditional variance model is included). If the number of rows 
%          exceeds the number necessary, then only the most recent observations
%          are used. If E0 is unspecified, any necessary presample observations 
%          are inferred from the corresponding presample response data Y0 (see
%          above), provided Y0 has at least Mdl.P + Mdl.Q observations, and 
%          presample predictor data X0 (if Mdl includes a regression component
%          and X0 is provided, see below); if Y0 is of insufficient length, 
%          presample observations of E0 are set to zero. The last row contains 
%          the most recent observation.
%
%   'V0'   Positive presample conditional variances, providing initial values 
%          for any conditional variance model; if the variance of the model 
%          is constant, then V0 is unnecessary. V0 is a column vector or a 
%          matrix. If V0 is a column vector, then it is applied to each 
%          forecasted path. If V0 is a matrix, then it must have numPaths 
%          columns (see notes below). V0 may have any number of rows, provided 
%          sufficient observations exist to initialize the variance model. 
%          If the number of rows exceeds the minimum, then only the most recent 
%          observations are used. If V0 is unspecified, any necessary presample 
%          observations are inferred from the corresponding residuals E0, 
%          provided E0 has sufficient observations required by the conditional 
%          variance model; if E0 is unspecified or of insufficient length, 
%          presample observations of V0 are set to the unconditional variance 
%          of the variance process. The last row contains the most recent 
%          observation.
%
%   'X0'   Presample matrix of predictor data used to infer presample 
%          innovations (see E0 above) when E0 is unspecified. Each column of 
%          X0 is a separate time series, and the last row of each contains 
%          the most recent observation of each series. To infer innovations 
%          in the absence of user-specified presample data E0, the number of 
%          observations in X0 must equal or exceed the number of observations 
%          in Y0 minus Mdl.P. When the number of observations exceeds the 
%          number necessary, only the most recent observations are used. 
%
%   'XF'   Matrix of forecasted (future) predictor data used to include a 
%          regression component in the model. XF represents the evolution of 
%          the predictor data in X0, forecasted into the future. The first 
%          row contains the 1-period-ahead forecast, the second row the 
%          2-period-ahead forecast, and so on. When the number of forecasts 
%          exceeds numPeriods, only the first numPeriods forecasts are used. 
%          If missing, the conditional mean will have no regression component 
%          regardless of the presence of any regression coefficients found 
%          in the model.
%
% Output Arguments:
%
%   Y - numPeriods-by-numPaths matrix of minimum mean square error (MMSE) 
%     forecasts of the response data Y0. The number of columns of Y (numPaths) 
%     is the largest number of columns of the presample arrays Y0, E0, and V0. 
%     The first row of Y contains the response forecasts in period 1, the 
%     second row contains the response forecasts in period 2, and so on until
%     the last row, which contains conditional mean forecasts at the forecast 
%     horizon. 
%
%   YMSE - numPeriods-by-numPaths matrix of mean square errors (MSE) of the
%     forecasts of the response Y. The number of columns of YMSE (numPaths) 
%     is the largest number of columns of the presample arrays Y0, E0, and 
%     V0. The first row of YMSE contains the forecast error variances in 
%     period 1, the second row contains the forecast error variances in 
%     period 2, and so on until the last row, which contains the forecast 
%     error variances at the forecast horizon. The square roots of YMSE are 
%     the standard errors of the forecasts Y above.
%
%   V - numPeriods-by-numPaths matrix of minimum mean square error (MMSE) 
%     forecasts of the conditional variances of future model innovations. The
%     number of columns of V (numPaths) is the largest number of columns of 
%     the presample arrays Y0, E0, and V0. The first row of V contains the 
%     conditional variance forecasts in period 1, the second row contains the
%     conditional variance forecasts in period 2, and so on until the last 
%     row, which contains conditional variance forecasts at the forecast 
%     horizon. 
%
% Notes:
%
%   o The number of sample paths (numPaths) is the largest column dimension 
%     of the presample arrays Y0, E0, and V0, but not fewer than one.
%
%   o If Y0, E0, and V0 are matrices with multiple columns (paths), they 
%     must have the same number of columns, otherwise an error occurs.
%
%   o Missing values, indicated by NaNs, are removed from Y0, E0, V0, and 
%     X0 by listwise deletion, thereby reducing the effective number of 
%     observations. That is, Y0, E0, V0, and X0 are merged into a composite 
%     series, and any row of the combined series with at least one NaN is 
%     removed. The presample data is assumed synchronized such that the last 
%     (most recent) observation of each series occurs at the same time. 
%     Likewise, any row of the forecasted predictor data XF with missing
%     observations is also removed.
%
%   o To include a regression component in the response forecast, the 
%     forecasted predictor data XF must be specified. Although the presample 
%     predictor data X0 is generally the same data referred to as 'X' used 
%     to estimate Mdl, X0 is used here only to infer presample innovations 
%     E0 when E0 is unspecified. When E0 is specified, X0 is ignored.
%     Although XF may be specified without X0, FORECAST issues an error when 
%     X0 is specified without XF.
%
%   o When computing the mean square errors of the conditional mean forecasts 
%     (see output YMSE above), the predictor data (see inputs X0 and XF 
%     above) are treated as exogenous, non-stochastic, and statistically 
%     independent of the model innovations. Therefore, the output YMSE 
%     reflects the variance associated with the ARIMA component of the input 
%     model Mdl alone.
%
% References:
%
%   [1] Baillie, R., and T. Bollerslev. "Prediction in Dynamic Models with 
%       Time-Dependent Conditional Variances." Journal of Econometrics.
%       Vol. 52, 1992, pp. 91-113.
%
%   [2] Bollerslev, T. "Generalized Autoregressive Conditional 
%       Heteroskedasticity." Journal of Econometrics. Vol. 31, 1986, pp.
%       307-327.
%
%   [3] Bollerslev, T. "A Conditionally Heteroskedastic Time Series Model
%       for Speculative Prices and Rates of Return." The Review Economics
%       and Statistics. Vol. 69, 1987, pp 542-547.
%
%   [4] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
%       Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
%       NJ: Prentice-Hall, 1994.
%
%   [5] Enders, W. Applied Econometric Time Series. Hoboken, NJ: John Wiley
%       & Sons, 1995.
%
%   [6] Engle, R. F. "Autoregressive Conditional Heteroskedasticity with
%       Estimates of the Variance of United Kingdom Inflation." 
%       Econometrica. Vol. 50, 1982, pp. 987-1007.
%
%   [7] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
% See also ARIMA, ESTIMATE, INFER, SIMULATE.

% Copyright 2018 The MathWorks, Inc.   

%
% Check input parameters and set defaults.
%

if numel(Mdl) > 1
   error(message('econ:arima:forecast:NonScalarModel'))
end

if nargin < 2
   error(message('econ:arima:forecast:NonEnoughInputs'))
end

%
% Determine if the user is calling the legacy interface in which 'Y0' is an
% optional N-V pair, or the new interface in which Y0 is a required input
% consistent with the newer VARM/FORECAST method.
%

X0     = zeros(0,size(Mdl.Beta,2));  % Default presample predictor data with correct # of columns
parser = inputParser;
parser.addRequired ('numPeriods',    @(x) validateattributes(x, {'double'}, {'scalar' 'integer' '>' 0}, '', 'forecast horizon'));
parser.addParameter('E0'        , 0, @(x) validateattributes(x, {'double'}, {}, '', 'presample residuals'));
parser.addParameter('V0'        , 0, @(x) validateattributes(x, {'double'}, {}, '', 'presample variances'));
parser.addParameter('X0'        ,X0, @(x) validateattributes(x, {'double'}, {}, '', 'presample regression matrix'));
parser.addParameter('XF'        ,[], @(x) validateattributes(x, {'double'}, {}, '', 'forecast regression matrix'));

try 
  if (nargin > 2) && isa(varargin{1}, 'double')  % Is it the new syntax?
     parser.addRequired ('Y0',    @(x) validateattributes(x, {'double'}, {}, '', 'presample responses'));
  else                                           % It's the old N-V pair syntax.
     parser.addParameter('Y0', 0, @(x) validateattributes(x, {'double'}, {}, '', 'presample responses'));
  end
  parser.parse(numPeriods, varargin{:});
catch exception
  exception.throwAsCaller();
end

horizon = parser.Results.numPeriods;
Y0      = parser.Results.Y0;
E0      = parser.Results.E0;
V0      = parser.Results.V0;
X0      = parser.Results.X0;
XF      = parser.Results.XF;

%
% Get model parameters and extract lags associated with non-zero coefficients.
%

constant = Mdl.PrivateConstant;                % Additive constant
variance = Mdl.PrivateVariance;                % Conditional variance

if any(strcmp(class(variance), {'garch' 'gjr' 'egarch'})) % Allow for a conditional variance model
   P = max(Mdl.PrivateVariance.P, Mdl.P);
   Q = max(Mdl.PrivateVariance.Q, Mdl.Q);                 % Total number of lagged e(t) needed
else
   P = Mdl.P;                                             % Total number of lagged y(t) needed
   Q = Mdl.Q;                                             % Total number of lagged e(t) needed
end

AR         = getLagOp(Mdl, 'Compound AR'); 
isARstable = isStable(AR);                      % Determine if the process is AR stable

AR     = reflect(AR);                           % This negates the AR coefficients
MA     = getLagOp(Mdl, 'Compound MA'); 
LagsAR = AR.Lags;                               % Lags of non-zero AR coefficients
LagsMA = MA.Lags;                               % Lags of non-zero MA coefficients
LagsAR = LagsAR(LagsAR > 0);                    % Exclude lag zero
LagsMA = LagsMA(LagsMA > 0);                    % Exclude lag zero

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
   error(message('econ:arima:forecast:UnspecifiedConstant'))
end

if any(isnan(AR))
   error(message('econ:arima:forecast:UnspecifiedAR'))
end

if any(isnan(MA))
   error(message('econ:arima:forecast:UnspecifiedMA'))
end

isX0specified = ~any(strcmpi('X0', parser.UsingDefaults));
isXFspecified = ~any(strcmpi('XF', parser.UsingDefaults));

if isX0specified && ~isXFspecified
   error(message('econ:arima:forecast:InconsistentRegressionComponent'))
end

if isXFspecified
   beta = Mdl.PrivateBeta;
   if isempty(beta) || any(isnan(beta))
      error(message('econ:arima:forecast:UnspecifiedBeta'))
   end
   if (numel(beta) ~= size(X0,2)) || (numel(beta) ~= size(XF,2))
      error(message('econ:arima:forecast:InconsistentRegression'))
   end
end

if ~any(strcmp(class(variance), {'garch' 'gjr' 'egarch'})) && any(isnan(variance))
   error(message('econ:arima:forecast:UnspecifiedVariance'))
end

%
% Compute the total number of observations generated for each path as the
% sum of the number of observations forested and the number of presample 
% observations needed to initialize the recursions.
%

maxPQ = max([P Q]);              % Maximum presample lags needed
T     = horizon + maxPQ;         % Total number of periods required for forecasting

%
% Compute the number of sample paths as the largest column dimension of the 
% presample arrays, but not fewer than one. If any of Y0, E0, and V0 are 
% column vectors, then the function "checkPresampleData" called later on will
% automatically expand the number of columns to the correct number of paths.
%

numPaths = max([size(Y0,2) size(E0,2) size(V0,2) 1]);

[nRows, nColumns] = size(Y0);

if nColumns ~= numPaths
   if (nColumns ~= 1) && ( ((nRows == 0) && (nColumns > 0)) || ~isempty(Y0) )
      error(message('econ:arima:forecast:InvalidY0', numPaths))
   end
end

[nRows, nColumns] = size(E0);

if nColumns ~= numPaths
   if (nColumns ~= 1) && ( ((nRows == 0) && (nColumns > 0)) || ~isempty(E0) )
      error(message('econ:arima:forecast:InvalidE0', numPaths))
   end
end

[nRows, nColumns] = size(V0);

if nColumns ~= numPaths
   if (nColumns ~= 1) && ( ((nRows == 0) && (nColumns > 0)) || ~isempty(V0) )
      error(message('econ:arima:forecast:InvalidV0', numPaths))
   end
end

%
% Remove missing observations (NaN's) via listwise deletion.
%

if any(isnan(Y0(:))) || any(isnan(X0(:))) || any(isnan(E0(:))) || any(isnan(V0(:)))
   [Y0, X0, E0, V0] = internal.econ.LagIndexableTimeSeries.listwiseDelete(Y0, X0, E0, V0);  % Pre-sample data
end

if any(isnan(XF(:)))
   XF = internal.econ.LagIndexableTimeSeries.listwiseDelete(XF);  % Forecast of predictor data
end

%
% Check any user-specified presample observations used for conditioning, or 
% generate any required observations automatically.
%

isY0specified = ~any(strcmpi('Y0', parser.UsingDefaults));
isE0specified = ~any(strcmpi('E0', parser.UsingDefaults));

if isE0specified      % Did the user specify presample e(t) observations?

%
%  Check user-specified presample data for the residuals e(t). 
%
%  Notice that the following line of code saves the original E0 input to 
%  forecast the conditional variance model later (rather than overwriting
%  it with a stripped version of itself).
%

   e0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(maxPQ,numPaths), 'E0', E0, Q);

%
%  Prepend the residuals with any user-specified presample observations and
%  transpose for efficiency.
%

   E            = zeros(numPaths,T);
   E(:,1:maxPQ) = e0';

else

%
%  The user did not specify presample e(t) observations. 
%

   if isY0specified && ( (size(Y0,1) >= (P + Q)) && ~isempty(Y0) )

      isE0Inferred = true;           % Initialize to true

      if isX0specified               % ARMAX model
%
%        Since the model has a regression component, whether or not sufficient 
%        observations of the input series y(t) required to infer initial values
%        of the residuals e(t) exist depends upon whether sufficient observations 
%        of the predictor data exist.
%
         if size(X0,1) >= (size(Y0,1) + Mdl.P)
%
%           When the number of presample observations of X0 equals or exceeds
%           the number of observations in Y0 plus Mdl.P, the INFER method is
%           called without specifying optional presample responses (see the 
%           Y0 input to INFER). When called in this manner, INFER automatically 
%           generates any required initial responses by backcasting, which
%           is consistent with the behavior in the absence of a regression
%           component (see just below).
%
            residuals = infer(Mdl, Y0, 'X', X0);  % Use backcasting

         else

            residuals = zeros(size(Y0));          % Pre-allocate to correct size
%
%           When the number of presample observations of X0 is less than the
%           number of observations in Y0 plus Mdl.P, the INFER method must 
%           be called by specifying optional presample responses (see the
%           Y0 input to INFER). To do this, the optional presample responses 
%           passed to INFER are obtained by stripping the first Mdl.P 
%           observations of the presample historical input series (the Y0 
%           input to this FORECAST method), thereby reducing the effective 
%           sample size. 
%
%           When presample stripping occurs, we must again test to ensure
%           that sufficient observations remain.
%
            if size(Y0((Mdl.P + 1):end,:),1) >= (P + Q)
%
%              Sufficient presample observations of y(t) exist to infer
%              presample innovations e(t), so initialize the first few
%              observations to the unconditional standard deviation.
%
               if any(strcmp(class(variance), {'garch' 'gjr' 'egarch'}))
                  residuals(1:Mdl.P,:) = repmat(sqrt(Mdl.PrivateVariance.UnconditionalVariance), Mdl.P, numPaths);
               else
                  residuals(1:Mdl.P,:) = repmat(sqrt(Mdl.PrivateVariance), Mdl.P, numPaths);
               end

               residuals((Mdl.P + 1):end,:) = infer(Mdl, Y0((Mdl.P + 1):end,:), 'Y0', Y0(1:Mdl.P,:), 'X', X0);
            else
               isE0Inferred = false;                % Insufficient observations
            end

         end

      else            % ARIMA model only
%
%        Since the model has no regression component, sufficient observations 
%        of the input series y(t) have been specified and so initial values 
%        of the innovations e(t) may be inferred.
%
         residuals = infer(Mdl, Y0);

      end

      E            = zeros(numPaths,T);
      E(:,1:maxPQ) = residuals((end - maxPQ + 1):end,:)';

   else
%
%    Insufficient observations of the input series y(t) have been specified,
%    so initialize any required presample observations with the unconditional 
%    mean of zero.

     E            = zeros(numPaths,T);   % Unconditional mean of e(t)
     isE0Inferred = false;

   end

end


if isY0specified  % Did the user specify presample y(t)?

%
%  Check user-specified presample data for the residuals e(t).
%

   Y0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(maxPQ,numPaths), 'Y0', Y0, Mdl.P);

%
%  Size the responses y(t) and initialize with specified data.
%

   Y = [Y0'  zeros(numPaths,horizon)];

else

%
%  The user did not specify presample y(t) observations. 
%

   if isARstable && (sum(AR) ~= 1) && ~isXFspecified
%
%     The model is AR-stable and without a regression component, so compute 
%     the unconditional (i.e., long-run) mean of the y(t) process directly 
%     from the parameters of the model and use it to initialize any required 
%     presample observations.
%
      average = constant / (1 - sum(AR));
      Y       = repmat([average(ones(1,maxPQ)) zeros(1,horizon)], numPaths, 1);

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
% Check any user-specified forecasted regression data for sufficient observations.
%

if isXFspecified
   try
%
%    The following code segment ensures we strip data from the end of the
%    X data forecast (i.e., we strip the most distant forecasts).
%
     XF = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(horizon, size(XF,2)), 'XF', XF(end:-1:1,:), horizon);
     XF = XF(end:-1:1,:);
   catch exception
     error(message('econ:arima:forecast:InsufficientXFRows'))
   end
end

%
% Apply iterative expectations one forecast step at a time. Such forecasts
% require that the process e(t) is a serially uncorrelated, zero-mean process 
% with a symmetric conditional probability distribution (see [1], pp. 94-95).
%

coefficients = [constant  AR  MA]';
I            = ones(numPaths,1);

if isXFspecified                     % ARIMAX model

   X = [zeros(maxPQ,size(XF,2)) ; XF].';

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

if nargout > 1     % Compute additional outputs only if necessary
%
%  Forecast the conditional variances.
%
   if any(strcmp(class(variance), {'garch' 'gjr' 'egarch'}))   % Conditional variance model

      isV0specified = ~any(strcmpi('V0', parser.UsingDefaults));

      if isE0specified && isV0specified

         V = forecast(variance, horizon, 'Y0', E0, 'V0', V0);

      elseif isE0specified

         V = forecast(variance, horizon, 'Y0', E0);

      elseif isV0specified

         if isE0Inferred
            V = forecast(variance, horizon, 'V0', V0, 'Y0', residuals);
         else
            V = forecast(variance, horizon, 'V0', V0);
         end

      else

         if isE0Inferred
            V = forecast(variance, horizon, 'Y0', residuals);
         else
            V = forecast(variance, horizon);
         end

     end

   else                              % Constant variance model

     V = variance(ones(horizon,numPaths));

   end

%
%  Compute variances of forecast errors of y(t) by converting the ARIMA 
%  model to its truncated infinite-degree MA representation.
%

   wState  = warning;                           % Save warning state
   cleanUp = onCleanup(@() warning(wState));    % Restore warning state
   
   warning('off', 'econ:LagOp:mldivide:WindowNotOpen')   
   warning('off', 'econ:LagOp:mldivide:WindowIncomplete')


   MA   = mldivide(getLagOp(Mdl, 'Compound AR'), ...
                   getLagOp(Mdl, 'Compound MA'), ...
                  'Degree', horizon - 1, 'RelTol', 0, 'AbsTol', 0);

   MA   = cell2mat(toCellArray(MA));
   MA   = [MA  zeros(1, horizon - numel(MA))];
   YMSE = toeplitz(MA.^2, [1 zeros(1, horizon - 1)]) * V;

end

%
% Remove the first max(P,Q) values used to initialize the variance forecast 
% such that the t-th row of V(t) is the t-period-ahead forecast of the 
% conditional variance, and transpose to a conventional time series format.
%

Y = Y(:,(maxPQ + 1):end)';

end
