function varargout = estimate(Mdl, YData, varargin)
%ESTIMATE Estimate ARIMA model parameters
%
% Syntax:
%
%   [EstMdl,EstParamCov,logL,info] = estimate(Mdl,Y)
%   [EstMdl,EstParamCov,logL,info] = estimate(Mdl,Y,param1,val1,...)
%
% Description:
%
%   Given an observed univariate time series, estimate the parameters of an
%   ARIMA model. The estimation process infers the residuals of the 
%   underlying response series and then fits the model to the response data 
%   via maximum likelihood.
%
% Input Arguments:
%
%   Mdl - ARIMA model specification object, as produced by the ARIMA 
%     constructor or ARIMA/ESTIMATE method.
%
%   Y - Response data whose residuals and conditional variances are inferred
%     and to which the model Mdl is fit. Y is a column vector, and therefore
%     a single path of the underlying series. The last observation of Y is 
%     the most recent.
%
% Optional Input Parameter Name/Value Pairs:
%
%   'Y0'         Presample response data, providing initial values for the 
%                model. Y0 is a column vector, and may have any number of 
%                rows, provided at least Mdl.P observations exist to 
%                initialize the model. If the number of rows exceeds Mdl.P, 
%                then only the most recent Mdl.P observations are used. If
%                Y0 is unspecified, any necessary observations are backcasted
%                (i.e., backward forecasted). The last row contains the most
%                recent observation.
%
%   'E0'         Mean-zero presample innovations, providing initial values 
%                for the model. E0 is a column vector, and may have any 
%                number of rows, provided sufficient observations exist to 
%                initialize the ARIMA model as well as any conditional 
%                variance model (the number of observations required is at 
%                least Mdl.Q, but may be more if a conditional variance 
%                model is included). If the number of rows exceeds the 
%                number necessary, then only the most recent observations
%                are used. If E0 is unspecified, any necessary observations 
%                are set to zero. The last row contains the most recent 
%                observation.
%
%   'V0'         Positive presample conditional variances, providing initial
%                values for any conditional variance model; if the variance 
%                of the model is constant, then V0 is unnecessary. V0 is a 
%                column vector, and may have any number of rows, provided 
%                sufficient observations exist to initialize the variance 
%                model. If the number of rows exceeds the number necessary, 
%                then only the most recent observations are used. If V0 is 
%                unspecified, any necessary observations are set to the 
%                average squared value of the inferred residuals. The last 
%                row contains the most recent observation.
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
%   'Options'    Optimization options created with OPTIMOPTIONS. If specified,
%                default optimization parameters are replaced by those in 
%                options. The default is an OPTIMOPTIONS object designed for 
%                the optimization function FMINCON, with 'Algorithm' = 'sqp' 
%                and 'ConstraintTolerance' = 1e-7. See documentation for 
%                OPTIMOPTIONS and FMINCON for details. 
%
%   'Constant0'  Scalar initial estimate of the constant of the model. If 
%                missing, an initial estimate is derived from standard time 
%                series techniques.
%
%   'AR0'        Vector of initial estimates of non-seasonal autoregressive
%                coefficients. The number of coefficients in AR0 must equal
%                the number of non-zero coefficients associated with the AR 
%                polynomial (excluding lag zero). If missing, initial 
%                estimates are derived from standard time series techniques.
%
%   'SAR0'       Vector of initial estimates of seasonal autoregressive
%                coefficients. The number of coefficients in SAR0 must equal
%                the number of non-zero coefficients associated with the
%                SAR polynomial (excluding lag zero). If missing, initial 
%                estimates are derived from standard time series techniques.
%
%   'MA0'        Vector of initial estimates of non-seasonal moving average
%                coefficients. The number of coefficients in MA0 must equal
%                the number of non-zero coefficients associated with the
%                MA polynomial (excluding lag zero). If missing, initial 
%                estimates are derived from standard time series techniques.
%
%   'SMA0'       Vector of initial estimates of seasonal moving average
%                coefficients. The number of coefficients in SMA0 must equal
%                the number of non-zero coefficients associated with the 
%                SMA polynomial (excluding lag zero). If missing, initial 
%                estimates are derived from standard time series techniques.
%
%   'Beta0'      Vector of initial estimates of the regression coefficients.
%                The number of coefficients in Beta0 must equal the number 
%                of columns in the predictor data matrix X (see above). If 
%                missing, initial estimates are derived from standard time 
%                series techniques.
%
%   'DoF0'       Scalar initial estimate of the degrees-of-freedom parameter
%                (used for t distributions only, and must exceed 2). If 
%                missing, the initial estimate is 10.
%
%   'Variance0'  A positive scalar initial variance estimate associated with 
%                a constant-variance model, or a cell vector of parameter 
%                name-value pairs of initial estimates associated with a 
%                conditional variance model. As a cell vector, the parameter 
%                names must be valid coefficients recognized by the variance 
%                model. If missing, initial estimates are derived from 
%                standard time series techniques.
%
%   'Display'    String vector or cell vector of character vectors indicating
%                what information to display in the command window. Values are:
%  
%                VALUE           DISPLAY
%  
%                o 'off'         No display to the command window. 
%  
%                o 'params'      Display maximum likelihood parameter 
%                                estimates, standard errors, and t statistics.
%                                This is the default.
%
%                o 'iter'        Display iterative optimization information.
%
%                o 'diagnostics' Display optimization diagnostics.
%
%                o 'full'        Display 'params', 'iter', and 'diagnostics'.
%
% Output Arguments:
%
%   EstMdl - An updated ARIMA model specification object containing the 
%     parameter estimates.
%
%   EstParamCov - Variance-covariance matrix associated with model parameters 
%     known to the optimizer. The rows and columns associated with any 
%     parameters estimated by maximum likelihood contain the covariances of 
%     the estimation errors; the standard errors of the parameter estimates 
%     are the square root of the entries along the main diagonal. The rows 
%     and columns associated with any parameters held fixed as equality 
%     constraints contain zeros. The covariance matrix is computed by the
%     outer product of gradients (OPG) method.
%
%   logL - Optimized loglikelihood objective function value.
%
%   info - Data structure of summary information with the following fields:
%
%     exitflag    - Optimization exit flag (see FMINCON)
%     options     - Optimization options (see OPTIMOPTIONS)
%     X           - Vector of final parameter/coefficient estimates
%     X0          - Vector of initial parameter/coefficient estimates
%
% Notes:
%
%   o Unspecified initial coefficient estimates are indicated by NaNs, which
%     are derived from standard time series techniques.
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
%   o The parameters known to the optimizer and included in EstParamCov are 
%     ordered as follows:
%
%       - Constant
%       - Non-zero AR coefficients at positive lags
%       - Non-zero SAR coefficients at positive lags
%       - Non-zero MA coefficients at positive lags
%       - Non-zero SMA coefficients at positive lags
%       - Regression coefficients (models with regression components only)
%       - Variance parameters (scalar for constant-variance models, vector
%         of additional parameters otherwise)
%       - Degrees-of-freedom (t distributions only)
%
%   o When 'Display' is specified, it takes precedence over the 'Diagnostics' 
%     and 'Display' selections found in the optimization 'Options' input.
%     However, when 'Display' is unspecified, all selections related to the
%     display of optimization information found in 'Options' are honored.
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
%   [3] Greene, W. H. Econometric Analysis. Upper Saddle River, NJ: 
%       Prentice Hall, 3rd Edition, 1997.
%
%   [4] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
% See also ARIMA, FORECAST, INFER, SIMULATE.

% Copyright 2019 The MathWorks, Inc.

if numel(Mdl) > 1
   error(message('econ:arima:estimate:NonScalarModel'))
end

%
% Ensure the input series is a non-empty vector.
%

if (nargin < 2) || isempty(YData) || ~iscolumn(YData)
   error(message('econ:arima:estimate:InvalidSeries'))
end

%
% Get model parameters and extract lags associated with non-zero coefficients. 
%

AR      = getLagOp(Mdl,  'AR');
MA      = getLagOp(Mdl,  'MA');
SAR     = getLagOp(Mdl, 'SAR');
SMA     = getLagOp(Mdl, 'SMA');

LagsAR  = AR.Lags;                    % Lags of non-zero  AR coefficients
LagsMA  = MA.Lags;                    % Lags of non-zero  MA coefficients
LagsSAR = SAR.Lags;                   % Lags of non-zero SAR coefficients
LagsSMA = SMA.Lags;                   % Lags of non-zero SMA coefficients
LagsAR  = LagsAR(LagsAR > 0);
LagsMA  = LagsMA(LagsMA > 0);
LagsSAR = LagsSAR(LagsSAR > 0);
LagsSMA = LagsSMA(LagsSMA > 0);

%
% Check input parameters and set defaults.
%

default = optimoptions('fmincon');

parser  = inputParser;
parser.addRequired ('requiredY',                        @(x) validateattributes(x, {'double'}, {}, '', 'response data'));
parser.addParameter('Y0'       , 0                    , @(x) validateattributes(x, {'double'}, {}, '', 'presample data'));
parser.addParameter('E0'       , 0                    , @(x) validateattributes(x, {'double'}, {}, '', 'presample residuals'));
parser.addParameter('V0'       , 0                    , @(x) validateattributes(x, {'double'}, {}, '', 'presample variances'));
parser.addParameter('X'        , []                   , @(x) validateattributes(x, {'double'}, {}, '', 'regression matrix'));
parser.addParameter('Constant0', nan                  , @(x) validateattributes(x, {'double'}, {'scalar'}, '', 'initial constant'))
parser.addParameter('AR0'      , nan(1,numel(LagsAR)) , @(x) validateattributes(x, {'double' 'cell'}, {}, '', 'initial AR coefficients'))
parser.addParameter('MA0'      , nan(1,numel(LagsMA)) , @(x) validateattributes(x, {'double' 'cell'}, {}, '', 'initial MA coefficients'))
parser.addParameter('SAR0'     , nan(1,numel(LagsSAR)), @(x) validateattributes(x, {'double' 'cell'}, {}, '', 'initial SAR coefficients'))
parser.addParameter('SMA0'     , nan(1,numel(LagsSMA)), @(x) validateattributes(x, {'double' 'cell'}, {}, '', 'initial SMA coefficients'))
parser.addParameter('Beta0'    , []                   , @(x) validateattributes(x, {'double' 'cell'}, {}, '', 'initial regression coefficients'))
parser.addParameter('DoF0'     , nan                  , @(x) validateattributes(x, {'double'}, {'scalar'}, '', 'initial degrees-of-freedom'))
parser.addParameter('options'  , default              , @(x) validateattributes(x, {'struct' 'optim.options.Fmincon'}, {'scalar'}, '', 'options'));
parser.addParameter('Variance0', nan                  , @(x) validateattributes(x, {'double' 'cell'}, {}, '', 'initial variance model coefficients'))
parser.addParameter('display'  , 'params'             , @internal.econ.displayCheck);

try
  parser.parse(YData, varargin{:});
catch exception
  exception.throwAsCaller();
end

YData           = parser.Results.requiredY;
Y0              = parser.Results.Y0;
E0              = parser.Results.E0;
V0              = parser.Results.V0;
XData           = parser.Results.X;
beta0           = parser.Results.Beta0(:)';
Constant0       = parser.Results.Constant0;
AR0             = parser.Results.AR0(:)';
MA0             = parser.Results.MA0(:)';
SAR0            = parser.Results.SAR0(:)';
SMA0            = parser.Results.SMA0(:)';
DoF0            = parser.Results.DoF0;
Variance0       = parser.Results.Variance0;
options         = parser.Results.options;
display         = lower(parser.Results.display);

% 
% Check presample data.
%

if ~iscolumn(E0) && ~isempty(E0)
   error(message('econ:arima:estimate:NonColumnVectorE0'))
end

if ~iscolumn(Y0) && ~isempty(Y0)
   error(message('econ:arima:estimate:NonColumnVectorY0'))
end

if ~iscolumn(V0) && ~isempty(V0)
   error(message('econ:arima:estimate:NonColumnVectorV0'))
end

%
% Remove missing observations (NaN's) via listwise deletion.
%

if any(isnan(Y0(:))) || any(isnan(E0(:))) || any(isnan(V0(:)))
   [Y0, E0, V0] = internal.econ.LagIndexableTimeSeries.listwiseDelete(Y0, E0, V0);    % Pre-sample data
end

if any(isnan(YData(:))) || any(isnan(XData(:)))
   [YData,XData] = internal.econ.LagIndexableTimeSeries.listwiseDelete(YData, XData); % In-sample data
end

%
% Get some additional model parameters. 
%

isVarianceConstant = ~any(strcmp(class(Mdl.PrivateVariance), {'garch' 'gjr' 'egarch'}));   % Is it a constant variance model?
isDistributionT    = strcmpi(Mdl.PrivateDistribution.Name, 'T');

constant = Mdl.PrivateConstant;       % Additive constant
variance = Mdl.PrivateVariance;       % Conditional variance

if ~isVarianceConstant                       % Allow for a conditional variance model
   P = max(Mdl.PrivateVariance.P, Mdl.P);
   Q = max(Mdl.PrivateVariance.Q, Mdl.Q);    % Total number of lagged e(t) needed
else
   P = Mdl.P;                         % Total number of lagged y(t) needed
   Q = Mdl.Q;                         % Total number of lagged e(t) needed
end

maxPQ = max([P Q]);                   % Maximum presample lags needed
T     = maxPQ + numel(YData);         % Total number of observations (including presample)

%
% Check to ensure the input series y(t) is long enough. ARIMA models require 
% the larger of the following:
%
% (1) Since all data series are prepended with "maxPQ" observations, the 
%     objective function needs at least "maxPQ" observations for backcasting 
%     y(t) in the event the user has not specified presample y(t) data. This 
%     check is just a "mechanical convenience", and not dictated by theory.
%
% (2) To compute initial guesses for model coefficients, the auto-covariance
%     calculation needs at least P + Q + 1 observations of y(t). The check 
%     is dictated the theory, assuming we use the B&J auto-covariance approach.
%

if numel(YData) < max(Mdl.P + Mdl.Q + 1, maxPQ)
   error(message('econ:arima:estimate:NotEnoughData'))
end

% 
% Check AR, MA, SAR, and SMA initial guesses.
%

if iscell(AR0)
   AR0 = cell2mat(AR0);
end

if iscell(MA0)
   MA0 = cell2mat(MA0);
end

if iscell(SAR0)
   SAR0 = cell2mat(SAR0);
end

if iscell(SMA0)
   SMA0 = cell2mat(SMA0);
end

%
% Validate any initial estimates provided by the user.
%

if isVarianceConstant
   if iscell(Variance0) && isscalar(Variance0)
      Variance0 = cell2mat(Variance0);
   end

   if ~isscalar(Variance0) || ~isnumeric(Variance0) || (Variance0 <= 0) 
      error(message('econ:arima:estimate:InvalidVariance0'))
   end
end

try
  Constant0 = internal.econ.validateUnivariateEstimates('Constant0', Constant0, 1             , NaN);
  AR0       = internal.econ.validateUnivariateEstimates('AR0'      , AR0      , numel(LagsAR) , 1:numel(LagsAR));
  MA0       = internal.econ.validateUnivariateEstimates('MA0'      , MA0      , numel(LagsMA) , 1:numel(LagsMA));
  SAR0      = internal.econ.validateUnivariateEstimates('SAR0'     , SAR0     , numel(LagsSAR), 1:numel(LagsSAR));
  SMA0      = internal.econ.validateUnivariateEstimates('SMA0'     , SMA0     , numel(LagsSMA), 1:numel(LagsSMA));
  DoF0      = internal.econ.validateUnivariateEstimates('DoF0'     , DoF0     , 1             , NaN);
catch exception
  exception.throwAsCaller();
end

%
% Validate regression component of the conditional mean model.
%

beta = Mdl.PrivateBeta;

if ~any(strcmpi('X', parser.UsingDefaults))
   isRegressionIncluded = true;

   if isempty(beta)
      beta            = nan(1,size(XData,2));
      Mdl.PrivateBeta = beta;
   else
      if numel(beta) ~= size(XData,2)
         error(message('econ:arima:estimate:InconsistentBeta'))
      end
   end
   if isempty(beta0)
      beta0 = nan(1,size(XData,2));
   else
      if numel(beta0) ~= size(XData,2)
         error(message('econ:arima:estimate:InconsistentBeta0'))
      end
   end
else
   isRegressionIncluded = false;
   Mdl.PrivateBeta      = zeros(1,0);
   beta                 = [];
   beta0                = [];
end

%
% Initialize some model parameters.
%

if isempty(LagsAR)
   AR = [];
else
   AR = AR.Coefficients;            % Lag Indexed Array
   AR = -[AR{LagsAR}];              % Non-zero AR coefficients (vector)
end

if isempty(LagsSAR)
   SAR = [];
else
   SAR = SAR.Coefficients;          % Lag Indexed Array
   SAR = -[SAR{LagsSAR}];           % Non-zero SAR coefficients (vector)
end

if isempty(LagsMA)
   MA  = [];
else
   MA  = MA.Coefficients;           % Lag Indexed Array
   MA  = [MA{LagsMA}];              % Non-zero MA coefficients (vector)
end

if isempty(LagsSMA)
   SMA = [];
else
   SMA = SMA.Coefficients;          % Lag Indexed Array
   SMA = [SMA{LagsSMA}];            % Non-zero MA coefficients (vector)
end

if isDistributionT
   DoF = Mdl.PrivateDistribution.DoF;  % Degrees-of-freedom (t distributions only)
else
   DoF = inf;
end

%
% Set optimization options.
%

if ~any(strcmpi('options', parser.UsingDefaults))  % Has the user specified options?
%
%  Optimization options must be either a data structure or an object of
%  type "optim.options.Fmincon" created by OPTIMOPTIONS. The following 
%  code segment only updates a structure, which may not have been created 
%  with OPTIMSET, thereby offering the user some additional flexibility.
%
%  On the other hand, if the user specified an options object, then the
%  type and corresponding properties are much more tightly controlled. In
%  this case, there can be no "missing" parameters, and the notion of 
%  updating a default object seems unnecessary. Therefore, options objects 
%  are honored in their entirety.
%
   if isa(options, 'struct')        % Is it a structure?
%
%     The user has specified an options structure, likely, but not necessarily,
%     by calling OPTIMSET. The following code creates a default options 
%     structure, and then updates it with any non-empty/non-missing values 
%     found in the structure specified by the user.
%
      default = optimset('fmincon');
      default = optimset(default, 'Algorithm', 'sqp', 'TolCon'     , 1e-7, ...
                                  'Display'  , 'off', 'Diagnostics', 'off');
      options = optimset(default, options);
   end

else

   options = optimoptions(options, 'Algorithm', 'sqp', 'ConstraintTolerance', 1e-7, ...
                                   'Display'  , 'off', 'Diagnostics', 'off');
end

%
% Allow the input display flag to take precedence over the relevant
% information found in the OPTIONS structure/object.
%

if ~any(strcmpi('display', parser.UsingDefaults))  % Did the user specify display

   if any(ismember({'diagnostics' , 'full'}, display))
      options.Diagnostics = 'on';
   else
      options.Diagnostics = 'off';
   end

   if any(ismember({'iter' , 'full'}, display))
      options.Display = 'iter';
   else
      options.Display = 'off';
   end

end

%
% Set the constraint tolerance and the value to trap objective function errors.
%

tolerance = 2 * options.TolCon;
trapValue = 1e20;

%
% Create a logical SOLVE vector to identify equality constraints for 
% individual parameters. These equalities are associated with non-NaN 
% coefficients found in the input model (Mdl). Elements of SOLVE = TRUE 
% indicate the corresponding coefficient found in the model is a NaN and must 
% be estimated, while elements of SOLVE = FALSE indicate the corresponding 
% model coefficient is not a NaN and is held fixed.
%

nCoefficients = 1 + numel(LagsAR) + numel(LagsSAR) + ...
                    numel(LagsMA) + numel(LagsSMA) + size(XData,2);

iAR   = 2:(2 + numel(LagsAR) - 1);
iSAR  = (numel(iAR) + 2):(numel(iAR) + numel(LagsSAR) + 1);
iMA   = (numel(iAR) + numel(iSAR) + 2):(numel(iAR) + numel(iSAR) + numel(LagsMA) + 1);
iSMA  = (numel(iAR) + numel(iSAR) + numel(iMA) + 2):(numel(iAR) + numel(iSAR) + numel(LagsMA) + numel(LagsSMA) + 1);
iBeta = (numel(iAR) + numel(iSAR) + numel(iMA) + numel(iSMA) + 2):nCoefficients;

if ~isVarianceConstant      % Allow for a conditional variance model

   [~,~,~,info] = estimate(variance, YData, 'Initialize', 'options', options);

   nTotal = nCoefficients + numel(info.solve);
   solve  = false(nTotal, 1);

   if isDistributionT
      solve((nCoefficients + 1):(nTotal - 1)) = info.solve(1:end-1);
      solve(nTotal)                           = info.solve(end);
   else
      solve((nCoefficients + 1):end) = info.solve;
   end

else                        % Constant variance model

   nTotal = nCoefficients + 1 + isDistributionT;
   solve  = false(nTotal, 1);

   if isDistributionT
      solve(nTotal - 1) = isnan(variance);
      solve(nTotal)     = isnan(Mdl.PrivateDistribution.DoF);
   else
      solve(nTotal) = isnan(variance);
   end

end

solve(1)     = isnan(constant);
solve(iAR)   = isnan(AR);
solve(iSAR)  = isnan(SAR);
solve(iMA)   = isnan(MA);
solve(iSMA)  = isnan(SMA);
solve(iBeta) = isnan(beta);

%
% Initialize the vector of initial estimates (X0) and the vector of 
% coefficients included in the model (X). Any element of X0 = NaN indicates 
% that the user did not specify an initial estimate for the corresponding
% coefficient, and will be updated later. Any element of X = NaN indicates
% that the corresponding coefficient will be estimated, and the corresponding
% element of SOLVE = TRUE; any non-NaN element of X is used as an initial 
% guess and copied to the corresponding element of X0, and the corresponding
% element of SOLVE = FALSE.
%

if ~isVarianceConstant 
   X  = [constant    AR    SAR    MA   SMA   beta   info.X']';
   X0 = [Constant0   AR0   SAR0   MA0  SMA0  beta0  info.X0']';
else
   X  = [constant    AR    SAR    MA   SMA   beta   variance   DoF(isDistributionT)]';
   X0 = [Constant0   AR0   SAR0   MA0  SMA0  beta0  Variance0  DoF0(isDistributionT)]';
end

X0(~solve) = X(~solve);  % Copy model coefficients to initial guesses

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
        XData = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(T,size(XData,2)), 'X', XData, numel(YData));
     else
%
%       No presample responses (Y0) are specified, and so the number of 
%       observations in X must equal or exceed the number of observations 
%       in Y plus the degree of the autoregressive polynomial (Mdl.P).
%
        XData = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(T,size(XData,2)), 'X', XData, numel(YData) + Mdl.P);
     end

   catch exception
     error(message('econ:arima:estimate:InsufficientXRows', numel(YData) + (Mdl.P * (~isY0specified))))
   end

end

%
% Check any user-specified presample observations used for conditioning, or 
% generate any required observations automatically.
%

if ~any(strcmpi('E0', parser.UsingDefaults))  % Did the user specify presample e(t) observations?

%
%  Check user-specified presample data for the residuals e(t).
%
   try
     E0              = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(maxPQ,1), 'E0', E0, Q);
     userSpecifiedE0 = true;
   catch exception
     exception.throwAsCaller();
   end

else

%
%  Set presample residuals e(t) = 0.
%

   E0              = zeros(maxPQ,1);
   userSpecifiedE0 = false;               % The user did not specify presample residuals e(t).

end

if isY0specified                          % Did the user specify presample y(t) observations?

%
%  Check user-specified presample data for the responses y(t).
%
   try
     Y0              = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(maxPQ,1), 'Y0', Y0, Mdl.P);
     userSpecifiedY0 = true;
   catch exception
     exception.throwAsCaller();
   end

else

%
%  Any presample observations are inferred during estimation, and so the 
%  presample vector must simply be pre-allocated correctly.
%

   Y0              = zeros(maxPQ,1);
   userSpecifiedY0 = false;               % The user did not specify presample y(t).

end

if ~any(strcmpi('V0', parser.UsingDefaults))  % Did the user specify presample v(t) observations?

%
%  Check user-specified presample data for the conditional variances v(t).
%
   try
     if isVarianceConstant
        V0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(ones(maxPQ,1), 'V0', V0, 0);
     else
        if isa(Mdl.PrivateVariance, 'egarch')
           V0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(ones(maxPQ,1), 'V0', V0, max(Mdl.PrivateVariance.P, Mdl.PrivateVariance.Q));
        else
           V0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(ones(maxPQ,1), 'V0', V0, Mdl.PrivateVariance.P);
        end
     end
     userSpecifiedV0 = true;
   catch exception
     exception.throwAsCaller();
   end

else

%
%  Any presample observations are inferred during estimation, and so the 
%  presample vector must simply be pre-allocated correctly.
%

   V0              = ones(maxPQ,1);
   userSpecifiedV0 = false;               % The user did not specify presample v(t).

end

%
% Provide initial estimates for coefficients associated with the ARIMA model
% only (i.e., including the variance of constant-variance models, but excluding 
% any conditional variance parameters and the degrees-of-freedom).
%

if any(solve(1:(nCoefficients + isVarianceConstant)))

%
%  Estimate parameters of an ARIMAX model by one of the following approaches:
%
%  (1) If the model has a regression component OR no MA component, then 
%      estimate the constant, AR, and predictor coefficients by OLS regression. 
%
%      After the ARX component is estimated, a subsequent estimation of any 
%      MA coefficients is made by applying the auto-covariance approach outlined 
%      in Box, Jenkins, and Reinsel (Appendix A6.2, pages 220-21), to the 
%      filtered residuals obtained in the first step, modified to account 
%      for any user-specified coefficient values.
%
%  (2) If the model has no regression component AND an MA component, then 
%      estimate the AR and MA coefficients by the auto-covariance method 
%      outlined in Box, Jenkins, and Reinsel (Appendix A6.2, pages 220-21), 
%      modified to account for any user-specified coefficient values.
%
   if isRegressionIncluded || (Mdl.Q == 0)
%
%     Create the non-seasonal and seasonal polynomials associated with
%     initial guesses of autoregressive coefficients.
%
      LagOpAR0  = LagOp([1 ; -X0(iAR) ], 'Lags', [0 LagsAR ]);
      LagOpSAR0 = LagOp([1 ; -X0(iSAR)], 'Lags', [0 LagsSAR]);
%
%     The lag operator polynomials just created exclude lags associated with 
%     zero-valued coefficient initial guesses (strictly speaking, coefficients 
%     whose magnitude <= LagOp.ZeroTolerance).
%
%     Therefore, the function ARX0 (called below) will then return only 
%     those non-zero autoregressive coefficients associated with positive 
%     lags found in the LagOpAR0 and LagOpSAR0 polynomials, with the AR 
%     coefficients followed by the SAR coefficients.
%
%     To properly account for a discrepancy in the event of zero-valued initial 
%     guesses, the following code segment pre-allocates the AR vector to the 
%     correct size, as indicated by the input model Mdl, with all zeros, then 
%     identifies which lags returned by the ARX0 function are to be inserted 
%     into the vector of initial guesses.
%
      AR       = zeros(numel(LagsAR) + numel(LagsSAR), 1);
      LagsAR0  = LagOpAR0.Lags;   
      LagsSAR0 = LagOpSAR0.Lags;
      LagsAR0  = LagsAR0(LagsAR0 > 0);        % Positive  AR lags of initial guesses
      LagsSAR0 = LagsSAR0(LagsSAR0 > 0);      % Positive SAR lags of initial guesses
      I1       = ismember(LagsAR , LagsAR0);  % AR  elements to fill
      I2       = ismember(LagsSAR, LagsSAR0); % SAR elements to fill
%
%     Create a modified model object to store the initial guesses, equality 
%     constraints, and NaNs indicating which coefficients to estimate.
%
      obj          = Mdl; 
      obj.Constant = X0(1);
      obj          = setLagOp(obj,  'AR', LagOpAR0 , 'Exclude Validation'); 
      obj          = setLagOp(obj, 'SAR', LagOpSAR0, 'Exclude Validation');
      obj.Beta     = X0(iBeta);
%
%     Estimate the ARX coefficients via OLS.
%
      [AR([I1 I2]),beta,constant,variance,residuals] = internal.econ.arx0(obj, YData, XData);

      MA = zeros(Mdl.Q, 1);

      if Mdl.Q > 0
%
%        Create the non-seasonal and seasonal product polynomial associated 
%        with initial guesses of moving average coefficients, and estimate 
%        via the auto-covariance approach of Box & Jenkins applied to the
%        filtered residuals derived from the preceding ARX estimation.
%
         LagOpMA0 = LagOp([1 ; X0(iMA) ], 'Lags', [0 LagsMA ]) * LagOp([1 ; X0(iSMA)], 'Lags', [0 LagsSMA]);
         [~, MA0] = internal.econ.arma0(residuals((Mdl.P + 1):end), LagOp(1), LagOpMA0);
%
%        Retain only those MA lags found in the model.
%
         LagsMA0     = LagOpMA0.Lags; 
         LagsMA0     = LagsMA0(LagsMA0 > 0);        % Positive MA lags of initial guesses
         MA(LagsMA0) = MA0(LagsMA0);                % Insert into the appropriate elements.
         MA          = [MA(LagsMA) ; MA(LagsSMA)];

      end

   else
%
%     Create the non-seasonal and seasonal product polynomials associated 
%     with initial guesses of autoregressive (excluding integration effects)
%     and moving average coefficients.
%
      LagOpAR0 = LagOp([1 ; -X0(iAR)], 'Lags', [0 LagsAR]) * LagOp([1 ; -X0(iSAR)], 'Lags', [0 LagsSAR]);
      LagOpMA0 = LagOp([1 ;  X0(iMA)], 'Lags', [0 LagsMA]) * LagOp([1 ;  X0(iSMA)], 'Lags', [0 LagsSMA]);
%
%     Create the product polynomial associated with non-seasonal and
%     seasonal integration, filter the input series y(t), then estimate.
%
      I = getLagOp(Mdl, 'Integrated Non-Seasonal') * getLagOp(Mdl, 'Integrated Seasonal');

      [AR0, MA0, constant, variance] = internal.econ.arma0(I(YData), LagOpAR0, LagOpMA0);
%
%     The function ARMA0, just called, returns all coefficients associated 
%     with the autoregressive and moving average product polynomials LagOpAR0
%     and LagOpMA0, respectively, at lags 1, 2, ... to the degrees of the
%     product polynomials.
%
%     However, the lag operator polynomials just created exclude lags 
%     associated with zero-valued coefficient initial guesses (strictly 
%     speaking, coefficients whose magnitude <= LagOp.ZeroTolerance), and 
%     must be accounted for.
%
%     The following code segment pre-allocates the resulting AR and MA product
%     vectors to the correct sizes, as indicated by the input model Mdl, with 
%     all zeros, then inserts any initial guesses into the appropriate elements.
%      
      AR          = zeros(Mdl.P - Mdl.PrivateD - Mdl.PrivateSeasonality, 1);
      MA          = zeros(Mdl.Q, 1);
      
      LagsAR0     = LagOpAR0.Lags;
      LagsMA0     = LagOpMA0.Lags; 
      LagsAR0     = LagsAR0(LagsAR0 > 0);        % Positive AR lags of initial guesses
      LagsMA0     = LagsMA0(LagsMA0 > 0);        % Positive MA lags of initial guesses

      AR(LagsAR0) = AR0(LagsAR0);                % Insert into the appropriate elements.
      MA(LagsMA0) = MA0(LagsMA0);                % Insert into the appropriate elements.
      AR          = [AR(LagsAR) ; AR(LagsSAR)];  % Retain only those lags found in the model
      MA          = [MA(LagsMA) ; MA(LagsSMA)];  % Retain only those lags found in the model
   
   end
%
%  Ensure the initial variance estimate is positive.
%
   if variance <= 0
      variance = tolerance;
   end
%
%  Any model coefficients to be estimated, and without any user-specified 
%  initial estimates, are overwritten with auto-generated initial estimates.
%
   i     = 1:(nCoefficients + isVarianceConstant);
   i     = isnan(X0(i)) & solve(i);
   x0    = [constant ; AR(:) ; MA(:) ; beta(:) ; variance(isVarianceConstant)]; 
   X0(i) = x0(i);

end

%
% Given the initial estimates of the ARIMA parameters, infer the residuals
% and estimate any parameters associated with a conditional variance model.
%

if ~isVarianceConstant      % Allow for a conditional variance model

   obj = setLagOp(Mdl, 'AR' , LagOp([1 -X0(iAR)' ], 'Lags', [0 LagsAR ]), 'Exclude Validation');
   obj = setLagOp(obj, 'SAR', LagOp([1 -X0(iSAR)'], 'Lags', [0 LagsSAR]), 'Exclude Validation');
   obj = setLagOp(obj, 'MA' , LagOp([1  X0(iMA)' ], 'Lags', [0 LagsMA ]), 'Exclude Validation');
   obj = setLagOp(obj, 'SMA', LagOp([1  X0(iSMA)'], 'Lags', [0 LagsSMA]), 'Exclude Validation');
%
%  Ensure the MA and SMA polynomials are invertible so the inferred residuals 
%  are not explosive.
%
   if ~isStable(getLagOp(obj,'MA'))
      obj = setLagOp(obj, 'MA' , LagOp([1  zeros(1,numel(iMA))] , 'Lags', [0 LagsMA]) , 'Exclude Validation');
   end
   if ~isStable(getLagOp(obj,'SMA'))
      obj = setLagOp(obj, 'SMA', LagOp([1  zeros(1,numel(iSMA))], 'Lags', [0 LagsSMA]), 'Exclude Validation');
   end

   obj.PrivateConstant = constant;
   obj.PrivateVariance = 1;               % Just a placeholder to avoid failure

   if isRegressionIncluded
      obj.PrivateBeta  = X0(iBeta)';
      residuals = infer(obj, YData, 'X', XData);
   else
      residuals = infer(obj, YData);
   end
   if iscell(Variance0)
      [~,~,~,info] = estimate(Mdl.PrivateVariance, residuals, 'Initialize', 'options', options, Variance0{:});
   else
      [~,~,~,info] = estimate(Mdl.PrivateVariance, residuals, 'Initialize', 'options', options);
   end

   X0((nCoefficients + 1):end) = info.X0;

end

%
% Identify the objective function and initialize common information.
%

T          = maxPQ + numel(YData); % Total number of observations (including presample)
E          = zeros(1,T);           % Preallocate residuals
E(1:maxPQ) = E0;                   % Initialize with presample data
V          = zeros(1,T);           % Preallocate conditional variances
V(1:maxPQ) = V0;                   % Initialize with presample data
YData      = [Y0 ; YData]';        % Augment with presample data and transpose 
XData      = XData';               % Transpose for consistency with other series

AR = getLagOp(Mdl, 'Compound AR'); % Compound AR polynomial needed for its non-zero lags
MA = getLagOp(Mdl, 'Compound MA'); % Compound MA polynomial needed for its non-zero lags

F  = @(X) nLogLike(X, YData, XData, E, V, Mdl, AR.Lags, MA.Lags, maxPQ, T, isDistributionT, options, ...
                   userSpecifiedY0, userSpecifiedE0, userSpecifiedV0, trapValue);
%
% Test to see if any coefficients are estimated. If no coefficients are
% estimated, then the entire estimation is a no-op!
%

if any(solve)  

%
% If necessary, update initial estimates for the degrees-of-freedom and set
% constraints for model coefficients.
%

  if isDistributionT && isnan(X0(nTotal))
     X0(nTotal) = 10;
  end

%
% Create the constraints for the ARIMA model and, if necessary, the
% constraints for any contained conditional variance model as well.
%
% In the following, the constraints and parameters known to FMINCON
% are ordered as follows: 
%
%    o Constant
%    o AR
%    o SAR
%    o MA
%    o SMA
%    o Beta regression coefficients (models with regression components only)
%    o Variance parameters (a scalar for constant-variance models, a vector
%      of additional parameters otherwise)
%    o Degrees-of-freedom (t distributions only)
%

  if isVarianceConstant      % Is it a constant-variance model?

%
%    Constant-variance ARIMA models handle the degrees-of-freedom parameter 
%    of t distributions directly.
%
     indices   = 1:(nCoefficients + 1 + isDistributionT);
     conStruct = internal.econ.arimaLinearConstraints(LagsAR, LagsSAR, LagsMA, LagsSMA, ...
                          beta, isDistributionT, isVarianceConstant, solve(indices), X0(indices), tolerance);

     if (numel(LagsAR) + numel(LagsSAR) + numel(LagsMA) + numel(LagsSMA)) == 0
         conStruct.nonLinear = [];  % No polynomials to remain stationary
     else
         conStruct.nonLinear = @(x) internal.econ.arimaNonLinearConstraints(x, LagsAR, LagsSAR, LagsMA, LagsSMA, tolerance);
     end

  else                       % Allow for a conditional variance model

%
%    Compute the linear constraints of the contained variance model and
%    augment the linear constraints of the ARIMA model with them.
%
%    ARIMA models with contained conditional variance models (e.g., GARCH) 
%    allow the contained variance model to handle the degrees-of-freedom 
%    parameter of t distributions.
%
     indices   = 1:nCoefficients;
     conStruct = internal.econ.arimaLinearConstraints(LagsAR, LagsSAR, LagsMA, LagsSMA, ...
                          beta, false, isVarianceConstant, solve(indices), X0(indices), tolerance);

     conStruct.lb  = [conStruct.lb ; info.constraints.lb];
     conStruct.ub  = [conStruct.ub ; info.constraints.ub];

     if isempty(info.constraints.A)
        conStruct.A = [];
        conStruct.b = [];
     else
        conStruct.A = [zeros(size(info.constraints.A,1),nCoefficients) info.constraints.A];
        conStruct.b = info.constraints.b;
     end
     if isempty(conStruct.Aeq)
        A_EQ1 = [];
     else
        A_EQ1 = [conStruct.Aeq  zeros(size(conStruct.Aeq,1),numel(info.X0))];
     end

     if isempty(info.constraints.Aeq)
        A_EQ2 = [];
     else
        A_EQ2 = [zeros(size(info.constraints.Aeq,1),nCoefficients)  info.constraints.Aeq];
     end

     conStruct.Aeq = [A_EQ1         ; A_EQ2];
     conStruct.beq = [conStruct.beq ; info.constraints.beq];

%
%    Compute the non-linear inequality constraints of the ARIMA model as well 
%    as the contained variance model if necessary.
%
     if isempty(info.constraints.nonLinear)
%
%       The contained variance model has no non-linear constraints, so just
%       use those of the ARIMA model if necessary.
%
        if (numel(LagsAR) + numel(LagsSAR) + numel(LagsMA) + numel(LagsSMA)) == 0
           conStruct.nonLinear = [];  % No polynomials to remain stable
        else
           conStruct.nonLinear = @(x) internal.econ.arimaNonLinearConstraints(x, LagsAR, LagsSAR, LagsMA, LagsSMA, tolerance);
        end

     else
%
%       The contained variance model has non-linear inequality constraints, 
%       so create a compound anonymous function to "concatenate" them. Also, 
%       wrap them up inside the DEAL function to ensure an empty matrix is
%       returned for the non-linear equality constraints, as required by
%       FMINCON.
%
        i = 1 + numel(LagsAR) + numel(LagsSAR) + numel(LagsMA) + numel(LagsSMA);
        f = @(x) [internal.econ.arimaNonLinearConstraints(x(1:i), LagsAR, LagsSAR, LagsMA, LagsSMA, tolerance) ; 
                                       info.constraints.nonLinear(x((i+1):end))];

        conStruct.nonLinear = @(x) deal(f(x),[]);  % FMINCON needs the 2nd output

     end

  end

%
% Estimate the coefficients.
%

  [coefficients, nLogL, exitflag] = fmincon(F, X0, conStruct.A, conStruct.b, ...
                 conStruct.Aeq, conStruct.beq, conStruct.lb, conStruct.ub, conStruct.nonLinear, options);

  if any(ismember({'params' , 'full'}, display))
     if any(coefficients <= conStruct.lb)
        warning(message('econ:arima:estimate:ActiveLowerBoundConstraints'));
     end
     if any(coefficients >= conStruct.ub)
        warning(message('econ:arima:estimate:ActiveUpperBoundConstraints'));
     end
     if ~isempty([conStruct.A conStruct.b]) && any(conStruct.A * coefficients >= conStruct.b)
        warning(message('econ:arima:estimate:ActiveLinearInequalityConstraints'));
     end
     if ~isempty(conStruct.nonLinear)
        if isVarianceConstant
           values = conStruct.nonLinear(coefficients);
        else
           if isempty(info.constraints.nonLinear)
              values = conStruct.nonLinear(coefficients);
           else
              [values,~] = conStruct.nonLinear(coefficients);
           end
        end
        if any(values >= 0)
           warning(message('econ:arima:estimate:ActiveNonLinearInequalityConstraints'));
        end
     end
  end

%
% It is possible that one or more of the coefficient estimates associated with 
% lag operator polynomials are sufficiently close to zero such that they 
% would be excluded from the output model. When this occurs, the output model 
% will not have the same non-zero lags as the input, and is usually an indication 
% that a simpler model is sufficient.
%
% Moreover, if the last coefficient associated with a lag operator polynomial
% is sufficiently close to zero then the degree P or Q of the input model will 
% differ from those of the output model. 
% 
% To ensure consistency between the degrees of the input and output models, the 
% following code segment identifies coefficients associated with lag operator 
% polynomials whose magnitudes are less than or equal to the LagOp zero-tolerance
% and sets them to LagOp.ZeroTolerance*2.
%
% Note for Future Enhancement: 
%
%    The value to which the magnitude of each polynomial coefficient is 
%    compared must be consistent with the 'Tolerance' parameter subsequently 
%    passed into any LagOp constructor below.
%

  if ~isempty(iAR) && (abs(coefficients(iAR(end))) <= LagOp.ZeroTolerance)
      coefficients(iAR(end)) = LagOp.ZeroTolerance*2;
  end

  if ~isempty(iSAR) && (abs(coefficients(iSAR(end))) <= LagOp.ZeroTolerance)
      coefficients(iSAR(end)) = LagOp.ZeroTolerance*2;
  end

  if ~isempty(iMA) && (abs(coefficients(iMA(end))) <= LagOp.ZeroTolerance)
      coefficients(iMA(end)) = LagOp.ZeroTolerance*2;
  end

  if ~isempty(iSMA) && (abs(coefficients(iSMA(end))) <= LagOp.ZeroTolerance)
      coefficients(iSMA(end)) = LagOp.ZeroTolerance*2;
  end

%
% Compute the variance-covariance matrix of the parameter estimates.
%

  try

    covariance = internal.econ.LagIndexableTimeSeries.errorCovarianceOpg(F, coefficients, solve, ...
                 conStruct.lb, conStruct.ub);

  catch exception
%
%   Indicate an error in the computation of the error covariance matrix 
%   by returning a matrix of NaN's.
%
    warning(message('econ:arima:estimate:CovarianceError'))
    covariance = nan(numel(solve));

  end
  
%
% Pack the estimated parameters into the output model.
%

  try
    Mdl = vectorToModel(Mdl, coefficients, covariance);
    
    Mdl.FitInformation.SampleSize             =  T - maxPQ;
    Mdl.FitInformation.LogLikelihood          = -nLogL;
    Mdl.FitInformation.NumEstimatedParameters =  sum(solve);
  catch old
    newMessage = message('econ:arima:estimate:InvalidVarianceModel');
    new        = MException(newMessage.Identifier,getString(newMessage));
    new        = addCause(new, old);
    throwAsCaller(new);
  end
%
% Pack relevant information into the output structure.
%
  information.exitflag = exitflag;
  information.options  = options;
  information.X        = coefficients;
  information.X0       = X0;

else

%
% The user entered a fully-specified model so just pack the private standard 
% errors information with zeros and the estimation indicators with FALSEs and 
% otherwise pass back the original model.
%
  covariance = zeros(numel(solve));
  nLogL      = F(X);
  Mdl        = vectorToModel(Mdl, X0, covariance);

  Mdl.FitInformation.SampleSize             =  T - maxPQ;
  Mdl.FitInformation.LogLikelihood          = -nLogL;
  Mdl.FitInformation.NumEstimatedParameters =  0;

  information.exitflag = [];
  information.options  = options;
  information.X        = X0;
  information.X0       = X0;

end

%
% Print estimation results.
%

if any(ismember({'full'}, display))
   summarize(Mdl)
elseif any(ismember({'params'}, display))
   disp(' ')
   disp("    " + getModelSummary(Mdl) + ":");
   disp(' ')
   results = summarize(Mdl); 
   disp(results.Table)
   if ~isVarianceConstant        % Conditional variance model
      disp(' ')
      disp(' ')
      disp("    " + getModelSummary(Mdl.PrivateVariance) + ":");
      disp(' ')
      disp(results.VarianceTable)
   end
end

varargout = {Mdl, covariance, -nLogL, information};

end  % Estimate method.

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function [nLogL,nLogLikelihoods] = nLogLike(X, YData, XData, E, V, Mdl, LagsAR, LagsMA, ...
          maxPQ, T, isDistributionT, options, isY0specified, isE0specified, isV0specified, trapValue)
%NLOGLIKE Negative log-likelihood function of univariate ARIMA processes
%
% Syntax:
%
%   [nLogL,nLogLikelihoods] = nLogLike(X, YData, XData, E, V, Mdl, LagsAR, ...
%               LagsMA, maxPQ, T, isDistributionT, options, isY0specified, ...
%               isE0specified, isV0specified, trapValue)
%
% Description:
%
%   Evaluate the negative log-likelihood objective function.
%
% Input Arguments:
%
%   X - Coefficient vector known to the optimizer.
%
%   YData - Pre-allocated row vector of observed data. The first maxPQ 
%     elements are initialized with any necessary pre-sample values provided 
%     by the user, or with zeros if no user-specified values are provided. 
%
%   XData - Matrix of predictor data used to determine if a regression 
%     component is included in the mean equation. If XData is empty, then 
%     the mean will have no regression component; if XData is not empty,
%     then it will be a pre-allocated, row-oriented data matrix with the same
%     number of columns as the response data in YData; the first maxPQ 
%     observations are initialized with zeros if user-specified response 
%     data are provided. 
%
%   E - Pre-allocated row vector of residuals. The first maxPQ elements
%     are initialized with any necessary pre-sample values provided by the 
%     user, or with zeros if no user-specified values are provided. 
%
%   V - Pre-allocated row vector of conditional variances. The first maxPQ 
%     elements are initialized with any necessary pre-sample values provided
%     by the user, or with zeros if no user-specified values are provided. 
%
%   Mdl - The ARIMA model whose coefficients are estimated.
%
%   LagsAR - Vector of non-zero lags of the compound AR polynomial of the 
%     ARIMA model.
%
%   LagsMA - Vector of non-zero lags of the compound MA polynomial of the 
%     ARIMA model.
%
%   maxPQ - The number of pre-sample observations padded to the beginning
%     of the input series YData, E, and V (see above).
%
%   T - The total number of observations of YData, E, and V.
%
%   isDistributionT - Logical flag: TRUE = if the distribution is 't', and
%     FALSE if the distribution is 'Gaussian'.
%
%   options - Optimization options structure (see OPTIMSET).
%
%   isY0specified - Logical flag: TRUE if the user specified pre-sample
%     observations of the time series YData, and FALSE otherwise. If FALSE,
%     then any necessary pre-sample observation of YData are backcasted.
%
%   isE0specified - Logical flag: TRUE if the user specified pre-sample
%     observations of the residuals E, and FALSE otherwise. If FALSE,
%     then any necessary pre-sample observation of E are set to zero for 
%     the ARIMA model but derived from the inferred residuals for any
%     conditional variance model contained in the 'Variance' property.
%
%   isV0specified - Logical flag: TRUE if the user specified pre-sample
%     observations of the conditional variances V, and FALSE otherwise. If 
%     FALSE, then any necessary pre-sample observation of V are set to zero 
%     for the ARIMA model but derived from the inferred residuals for any
%     conditional variance model contained in the 'Variance' property.
%
%   trapValue - The value assigned to the negative log-likelihood in the 
%     event its computed value is non-finite (i.e., Inf and NaN). This may 
%     occur in presence of non-finite residuals that arise due to unstable 
%     polynomials, or to excessively large log-likelihood components. This 
%     value is used to "trap" such conditions, and in the future may be 
%     dependent on the optimization algorithm (e.g., NaN for all algorithms
%     except 'Active-Set').
%
% Output Arguments:
%
%   nLogL - The value of the negative log-likelihood objective function.
%
%   nLogLikelihoods - Vector of individual terms which sum to the negative
%     log-likelihood (see nLogL above). The length of his vector will be 
%     T - maxPQ, and is needed to support the calculation of the parameter 
%     estimation error covariance matrix.
%

%
% Allow the predictor data to determine the presence of a regression component.
%

isRegressionIncluded = ~isempty(XData);

%
% Since the objective is to infer residuals, retain all the coefficients of 
% non-zero lags of the compound AR polynomial (including 0), but exclude the 
% zero-lag coefficient of the compound MA polynomial.
%

LagsMA = LagsMA(LagsMA > 0);            % Exclude zero-lag MA coefficient

%
% Extract component polynomials.
%

AR  = getLagOp(Mdl, 'AR');
MA  = getLagOp(Mdl, 'MA');
SAR = getLagOp(Mdl, 'SAR');
SMA = getLagOp(Mdl, 'SMA');

%
% Pack the polynomial coefficients known to the optimizer into the appropriate 
% lag of the corresponding component polynomial. As we pack the polynomials, 
% notice that zero-lag AR & SAR coefficients, defined to be 1, are excluded.
%

L_AR = AR.Lags;                         % Non-seasonal AR
L_AR = L_AR(L_AR > 0);                  % Exclude zero-lag coefficient

i1 = 2;
i2 = i1 + numel(L_AR) - 1;

AR = [1 zeros(1,AR.Degree)];
AR(L_AR + 1) = -X(i1:i2);

L_SAR = SAR.Lags;                       % Seasonal AR
L_SAR = L_SAR(L_SAR > 0);               % Exclude zero-lag coefficient

i1  = i2 + 1;
i2  = i1 + numel(L_SAR) - 1;

SAR = [1 zeros(1,SAR.Degree)];
SAR(L_SAR + 1) = -X(i1:i2);

L_MA = MA.Lags;                         % Non-seasonal MA
L_MA = L_MA(L_MA > 0);                  % Exclude zero-lag coefficient

i1  = i2 + 1;
i2  = i1 + numel(L_MA) - 1;

MA = [1 zeros(1,MA.Degree)];
MA(L_MA + 1) = X(i1:i2);

L_SMA = SMA.Lags;                       % Seasonal MA
L_SMA = L_SMA(L_SMA > 0);               % Exclude zero-lag coefficient

i1  = i2 + 1;
i2  = i1 + numel(L_SMA) - 1;

SMA = [1 zeros(1,SMA.Degree)];
SMA(L_SMA + 1) = X(i1:i2);

%
% Pack the regression coefficients into a dedicated coefficient vector.
%

if isRegressionIncluded
   i1   = i2 + 1;
   i2   = i1 + size(XData,1) - 1;
   beta = X(i1:i2)';
end

%
% Re-constitute the compound AR and MA polynomials and form the data 
% coefficient vector.
%

AR = conv(AR, SAR);

if Mdl.PrivateD > 0
   AR = conv(AR, cell2mat(toCellArray(getLagOp(Mdl, 'Integrated Non-Seasonal'))));
end

if Mdl.PrivateSeasonality > 0
   AR = conv(AR, cell2mat(toCellArray(getLagOp(Mdl, 'Integrated Seasonal'))));
end

MA = conv(MA, SMA);
MA = MA(2:end);                         % We solve for e(t), so exclude lag 0

%
% Initialize pre-sample y(t) data if necessary. 
%

if ~isY0specified && (Mdl.P > 0)

%
%  Backward-forecast (i.e., backcast) the observed data y(t) assuming all 
%  prior residuals e(t) = 0: first time-reverse the relevant segment of 
%  y(t), then forecast the time-reversed series, and finally time-reverse
%  the forecast.
%

   Y0           = [YData((2 * maxPQ):-1:(maxPQ + 1))  zeros(1,maxPQ)];
   coefficients = [X(1)  -AR(LagsAR(2:end) + 1)]';

   if isRegressionIncluded                     % ARMAX model
%
%     For fully-specified ARIMAX models, the following is equivalent to:
%
%     Y0 = forecast(Mdl, maxPQ, 'Y0', YData(end:-1:(maxPQ + 1))'  , ...
%                               'X0', XData(:,end:-1:(maxPQ + 1))', ...
%                               'XF', XData(:,maxPQ:-1:1)'        , 
%                               'E0', zeros(maxPQ,1));
%     YData(1:maxPQ) = Y0(end:-1:1);
%
      XData0 = XData(:,(2 * maxPQ):-1:1);

      for t = (maxPQ + 1):(2 * maxPQ)
          data  = [1  Y0(t - LagsAR(2:end))];  % Exclude lag 0 from the backcast
          Y0(t) = data * coefficients  +  beta * XData0(:,t);
      end

   else                                        % ARIMA model only
%
%     For fully-specified ARIMA models, the following is equivalent to:
%
%     Y0             = forecast(Mdl, maxPQ, 'Y0', YData(end:-1:(maxPQ + 1))', 'E0', zeros(maxPQ,1)); 
%     YData(1:maxPQ) = Y0(end:-1:1);
%
      for t = (maxPQ + 1):(2 * maxPQ)
          data  = [1  Y0(t - LagsAR(2:end))];  % Exclude lag 0 from the backcast
          Y0(t) = data * coefficients;
      end
   end

   YData(maxPQ:-1:1) = Y0((maxPQ + 1):end);

end

%
% Infer the residuals.
%

coefficients = [-X(1)  AR(LagsAR+1)  -MA(LagsMA)]';

if isRegressionIncluded                     % ARIMAX model

   E = internal.econ.arimaxMex(LagsAR, LagsMA, coefficients, maxPQ, YData', E', XData', -beta)'; 

%
%  Note:
%
%  To improve runtime performance, the following code segment has been replaced 
%  by the private/undocumented C-MEX file called above, which does exactly the 
%  same thing. To restore the original MATLAB code, simply uncomment the 
%  following code segment:


%    for t = (maxPQ + 1):T
%        data = [1  YData(t - LagsAR)  E(t - LagsMA)];
%        E(t) = data * coefficients  -  beta * XData(:,t);
%    end

else                                        % ARIMA model only

   E = internal.econ.arimaxMex(LagsAR, LagsMA, coefficients, maxPQ, YData', E')'; 

%
%  Note:
%
%  To improve runtime performance, the following code segment has been replaced 
%  by the private/undocumented C-MEX file called above, which does exactly the 
%  same thing. To restore the original MATLAB code, simply uncomment the 
%  following code segment:


%    for t = (maxPQ + 1):T
%        data = [1  YData(t - LagsAR)  E(t - LagsMA)];
%        E(t) = data * coefficients;
%    end
   
end

%
% Trap the presence of non-finite residuals (i.e., Inf's and NaN's) that may 
% arise due to unstable polynomials, and assign an appropriate log-likelihood value.
%

if any(~isfinite(E))

   nLogLikelihoods = nan(1, T - maxPQ);
   nLogL           = trapValue;

else

%
%  Evaluate the negative log-likelihood objective function.
%

   if any(strcmp(class(Mdl.PrivateVariance), {'garch' 'gjr' 'egarch'}))
%
%     The optimization involves a conditional variance model contained inside
%     the ARIMA model, so gather the required information without optimizing,
%     then evaluate the objective function.
%
      if isE0specified && isV0specified
         [~,~,~,info] = estimate(Mdl.PrivateVariance, E((maxPQ + 1):T)', 'Initialize', 'options', options, ...
                                'E0', E(1:maxPQ)', 'V0', V(1:maxPQ)');
      elseif isE0specified
         [~,~,~,info] = estimate(Mdl.PrivateVariance, E((maxPQ + 1):T)', 'Initialize', 'options', options, ...
                                'E0', E(1:maxPQ)');
         if (numel(info.E0) > 0) || (numel(info.V0) > 0)  % Then the following assignments will not fail 
            V(1:maxPQ) = info.V0(1);
         end
      elseif isV0specified
         [~,~,~,info] = estimate(Mdl.PrivateVariance, E((maxPQ + 1):T)', 'Initialize', 'options', options, ...
                                'V0', V(1:maxPQ)');
         if (numel(info.E0) > 0) || (numel(info.V0) > 0)  % Then the following assignments will not fail
            E(1:maxPQ) = info.E0(1);
         end
      else
         [~,~,~,info] = estimate(Mdl.PrivateVariance, E((maxPQ + 1):T)', 'Initialize', 'options', options);
         if (numel(info.E0) > 0) || (numel(info.V0) > 0)  % Then the following assignments will not fail
            E(1:maxPQ) = info.E0(1);
            V(1:maxPQ) = info.V0(1);
         end
      end

      if isDistributionT
         [nLogL,nLogLikelihoods] = Mdl.PrivateVariance.nLogLikeT(X((end - numel(info.X0) + 1):end), V, E, ...
                                       info.Lags, 1, maxPQ, T, info.presampleFlags, trapValue);
      else
         [nLogL,nLogLikelihoods] = Mdl.PrivateVariance.nLogLikeGaussian(X((end - numel(info.X0) + 1):end), V, E, ...
                                       info.Lags, 1, maxPQ, T, info.presampleFlags, trapValue);
      end

   else                         % Constant variance model

      if isDistributionT
%
%        Evaluate standardized Student's t negative log-likelihood function.
%
         V((maxPQ + 1):T) = X(end - 1);  % Assign the constant variance
         DoF              = X(end);      % Assign the degrees-of-freedom

         if DoF <= 200
% 
%           Standardized t-distributed residuals.
%
            nLogLikelihoods =  0.5 * (log(V((maxPQ + 1):T))  +  (DoF + 1) * log(1 + (E((maxPQ + 1):T).^2)./(V((maxPQ + 1):T) * (DoF - 2))));
            nLogLikelihoods =  nLogLikelihoods - log(gamma((DoF + 1)/2) / (gamma(DoF/2) * sqrt(pi * (DoF - 2))));
            nLogL           =  sum(nLogLikelihoods, 2);

         else
% 
%           Gaussian residuals.
%
            nLogLikelihoods = 0.5 * (log(2*pi*V((maxPQ + 1):T)) + ((E((maxPQ + 1):T).^2)./V((maxPQ + 1):T)));
            nLogL           = sum(nLogLikelihoods, 2);

         end

      else

         V((maxPQ + 1):T) = X(end);  % Assign the constant variance
%
%        Evaluate the Gaussian negative log-likelihood objective.
%
         nLogLikelihoods = 0.5 * (log(2*pi*V((maxPQ + 1):T)) + ((E((maxPQ + 1):T).^2)./V((maxPQ + 1):T)));
         nLogL           = sum(nLogLikelihoods, 2);

      end

%
%     Trap any degenerate negative log-likelihood values.
%

      i = isnan(nLogL) | (imag(nLogL) ~= 0) | isinf(nLogL);

      if any(i)
         nLogL(i) = trapValue;
      end

   end

end

end
