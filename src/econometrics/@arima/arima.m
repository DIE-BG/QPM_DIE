classdef (Sealed) arima < matlab.mixin.CustomDisplay
%ARIMA Create an ARIMA model
%
% Syntax:
%
%   Mdl = arima(p,D,q)
%   Mdl = arima(param1,val1,param2,val2,...)
%
% Description:
%
%   Create an ARIMA(p,D,q) model by specifying either the degrees p, D, and
%   q (short-hand syntax for non-seasonal models) or a list of parameter 
%   name-value pairs (long-hand syntax). For response process y(t) and model 
%   innovations e(t), either syntax allows the creation of an ARIMA model of 
%   the general form
%
%   y(t) = c + a1*y(t-1) + ... + aP*y(t-P) + e(t) + b1*e(t-1) + ... + bQ*e(t-Q)
%
%   with P auto-regressive terms and Q moving average terms (see notes below
%   for a discussion of the relationship between inputs and model degrees).
%
% Input Arguments (Short-Hand Syntax for Non-Seasonal Models Only):
%
%   p - Nonnegative integer indicating the degree of the non-seasonal 
%     auto-regressive polynomial.
%
%   D - Nonnegative integer indicating the degree of the non-seasonal
%     differencing polynomial (the degree of non-seasonal integration).
%
%   q - Nonnegative integer indicating the degree of the non-seasonal 
%     moving average polynomial.
%
% Input Arguments (Parameter Name/Value Pairs):
%
%   'Constant'  Scalar constant c of the model. If unspecified, the constant 
%               is set to NaN.
%
%   'AR'        A cell vector of non-seasonal auto-regressive coefficients.
%               When specified without corresponding lags, AR is a cell 
%               vector of coefficients at lags 1, 2, ... to the degree of 
%               the non-seasonal auto-regressive polynomial. When specified 
%               along with ARLags (see below), AR is a commensurate length 
%               cell vector of coefficients associated with the lags in 
%               ARLags. If unspecified, AR is a cell vector of NaNs the
%               same length as ARLags (see below).
%
%   'MA'        A cell vector of non-seasonal moving average coefficients.
%               When specified without corresponding lags, MA is a cell 
%               vector of coefficients at lags 1, 2, ... to the degree of 
%               the non-seasonal moving average polynomial. When specified 
%               along with MALags (see below), MA is a commensurate length 
%               cell vector of coefficients associated with the lags in
%               MALags. If unspecified, MA is a cell vector of NaNs the
%               same length as MALags (see below).
%
%   'ARLags'    A vector of positive integer lags associated with the AR
%               coefficients. If unspecified, ARLags is a vector of integers 
%               1, 2, ... to the degree of the non-seasonal auto-regressive 
%               polynomial (see AR above).
%
%   'MALags'    A vector of positive integer lags associated with the MA
%               coefficients. If unspecified, MALags is a vector of integers
%               1, 2, ... to the degree of the non-seasonal moving average 
%               polynomial (see MA above).
%
%   'SAR'       A cell vector of seasonal auto-regressive coefficients. When 
%               specified without corresponding lags, SAR is a cell vector 
%               of coefficients at lags 1, 2, ... to the degree of the 
%               seasonal auto-regressive polynomial. When specified along 
%               with SARLags (see below), SAR is a commensurate length cell 
%               vector of coefficients associated with the lags in SARLags. 
%               If unspecified, SAR is a cell vector of NaNs the same length 
%               as SARLags (see below).
%
%   'SMA'       A cell vector of seasonal moving average coefficients. When 
%               specified without corresponding lags, SMA is a cell vector 
%               of coefficients at lags 1, 2, ... to the degree of the 
%               seasonal moving average polynomial. When specified along 
%               with SMALags (see below), SMA is a commensurate length cell 
%               vector of coefficients associated with the lags in SMALags. 
%               If unspecified, SMA is a cell vector of NaNs the same length 
%               as SMALags (see below). 
%
%   'SARLags'   A vector of positive integer lags associated with the SAR
%               coefficients. If unspecified, SARLags is a vector of integers 
%               1, 2, ... to the degree of the seasonal auto-regressive 
%               polynomial (see SAR above)
%
%   'SMALags'   A vector of positive integer lags associated with the SMA
%               coefficients. If unspecified, SMALags is a vector of integers
%               1, 2, ... to the degree of the seasonal moving average
%               polynomial (see SMA above)
%
%   'D'         A nonnegative integer indicating the degree of the non-
%               seasonal differencing polynomial (the degree of non-seasonal 
%               integration). If unspecified, the default is zero.
%
%'Seasonality'  A nonnegative integer indicating the degree of the seasonal 
%               differencing polynomial. If unspecified, the default is
%               zero.
%
%   'Beta'      A vector of regression coefficients associated with a
%               regression component in the conditional mean. The presence 
%               of a non-empty Beta vector allows for an ARIMAX model in 
%               which predictor data is included in the ARIMA model equation
%               shown above.
%
%   'Variance'  A positive scalar variance of the model innovations, or a
%               supported conditional variance model object (e.g., a GARCH
%               model). If unspecified, the default is NaN.
%
%'Distribution' The conditional probability distribution of the innovations
%               process. Distribution is a string or character vector specified
%               as 'Gaussian' or 't', or a structure with field 'Name' which 
%               stores the distribution 'Gaussian' or 't'. If the distribution
%               is  't', then the structure must also have the field 'DoF' 
%               to store the degrees-of-freedom.
%
% 'Description' Description of the model, specified as a string. The default 
%               is a summary of the parametric form of the model.
%
% Output Argument:
%
%   Mdl - An ARIMA model with the following properties:
%
%         o Description
%         o Distribution
%         o P
%         o D
%         o Q
%         o Constant
%         o AR
%         o MA
%         o SAR
%         o SMA
%         o Beta
%         o Variance
%         o Seasonality
%
% Notes:
%
%   o The properties P and Q of ARIMA models do not necessarily conform to
%     standard Box and Jenkins notation. 
%
%     The property P is the degree of the compound autoregressive polynomial, 
%     or the total number of lagged observations of the underlying process 
%     necessary to initialize the auto-regressive component of the model. 
%     P includes the effects of non-seasonal and seasonal integration 
%     captured by the properties D and Seasonality, respectively, and the
%     non-seasonal and seasonal auto-regressive polynomials AR and SAR, 
%     respectively.
%
%     The property Q is the degree of the compound moving average polynomial, 
%     or the total number of lagged innovations of the underlying process 
%     necessary to initialize the moving average component of the model. Q 
%     includes the effects of non-seasonal and seasonal moving average 
%     polynomials MA and SMA, respectively.
%
%     If the model has no integration and no seasonal components, only then 
%     will the properties P and Q conform to standard Box and Jenkins 
%     notation for an ARIMA(P,0,Q) = ARMA(P,Q) model.
%
%   o The lags associated with the seasonal polynomials SAR and SMA are
%     specified in the periodicity of the observed data (as are AR and MA), 
%     and not as multiples of the Seasonality parameter. This convention
%     does not conform to standard Box and Jenkins notation, yet is a more
%     flexible approach to incorporate multiplicative seasonal models.
%
%   o The coefficients AR, SAR, MA, and SMA are each associated with an
%     underlying lag operator polynomial and subject to a near-zero 
%     tolerance exclusion test. That is, each coefficient is compared to 
%     the default zero tolerance 1e-12, and is included in the model only 
%     if the magnitude is greater than 1e-12; if the coefficient magnitude 
%     is less than or equal to 1e-12, then it is sufficiently close to zero 
%     and excluded from the model. See LagOp for additional details.
%
% See also ESTIMATE, FORECAST, INFER, SIMULATE, FILTER, IMPULSE.

% Copyright 2019 The MathWorks, Inc.

properties (GetAccess = public, SetAccess = private, Dependent)
  P                 % Degree of the composite autoregressive polynomial
  Q                 % Degree of the composite moving average polynomial
end

properties (GetAccess = public, SetAccess = private, Dependent, Hidden)
  UnconditionalMean % Unconditional (long-run) mean
end

properties (Access = public)
  Description = ''; % Model description/summary string
end

properties (Access = public, Dependent)
  Distribution     % Conditional probability distribution
  Constant         % Model constant
  AR               % Non-seasonal autoregressive coefficients
  SAR              % Seasonal autoregressive coefficients
  MA               % Non-seasonal moving average coefficients
  SMA              % Seasonal moving average coefficients
  D                % Degree of non-seasonal integration
  Seasonality      % Degree of seasonal integration
  Beta             % Regression coefficients
  Variance         % Model variance
end

properties (Access = private)
  PrivateDistribution = struct('Name', "Gaussian");
  PrivateConstant     = NaN;        % Model constant
  LagOpRHS            = {[] []};    % Lag operator polynomial of the RHS coefficients
  LagOpLHS            = {[] []};    % Lag operator polynomial of the LHS coefficients
  PrivateD            = 0;          % Degree of non-seasonal integration
  PrivateSeasonality  = 0;          % Degree of seasonal integration
  PrivateVariance     = NaN;        % Model variance
  PrivateBeta         = zeros(1,0); % Regression coefficients
  FitInformation      = [];         % Estimation information
end

methods    % GET/SET methods


%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function P = get.P(Mdl)
     P = Mdl.LagOpLHS{1}.Degree + Mdl.LagOpLHS{2}.Degree + Mdl.PrivateD + Mdl.PrivateSeasonality;
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function Q = get.Q(Mdl)
     Q = Mdl.LagOpRHS{1}.Degree + Mdl.LagOpRHS{2}.Degree;
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function value = get.Distribution(Mdl)
     value = Mdl.PrivateDistribution;
  end

  function Mdl = set.Distribution(Mdl, value)

   if ischar(value) || ( isstring(value) && isscalar(value) )
      if strcmpi(value, "Gaussian")
         value = struct('Name', "Gaussian");
      elseif strcmpi(value, "t")
         value = struct('Name', "t", 'DoF', NaN);
      else
         error(message('econ:arima:arima:InvalidDistributionString'))
      end

   elseif isstruct(value)
      fields = fieldnames(value);

      if ~any(strcmp('Name', fields))
         error(message('econ:arima:arima:MissingName'))
      end

      if strcmpi(value.Name, "Gaussian")
         if numel(fields) > 1
            error(message('econ:arima:arima:InvalidGaussian'))
         end
         value.Name = "Gaussian";

      elseif strcmpi(value.Name, "t")
         if ~any(strcmp('DoF', fields)) || (numel(fields) > 2)
            error(message('econ:arima:arima:InvalidT'))
         end

         if ~isscalar(value.DoF) || ~isa(value.DoF, 'double') || (value.DoF <= 2)
            error(message('econ:arima:arima:InvalidDoF'))
         end
         value.Name = "t";

      else
         error(message('econ:arima:arima:InvalidDistributionName'))
      end
      
   else
      error(message('econ:arima:arima:InvalidDistributionDataType'))
   end

   updateDescription       = Mdl.Description == getModelSummary(Mdl);
   Mdl.PrivateDistribution = value;
   Mdl.FitInformation      = [];           % Indicate that the model is NOT estimated.

   if ~isa(Mdl.PrivateVariance, 'double')  % Update any conditional variance model.
      Mdl.PrivateVariance.Distribution = value;
   end

   if updateDescription
      Mdl.Description = getModelSummary(Mdl);
   end

  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function Mdl = set.Description(Mdl, value)
     validateattributes(value, {'char'  'string'}, {'scalartext'}, '', 'Description');
     Mdl.Description = string(value);
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function value = get.Constant(Mdl)
     value = Mdl.PrivateConstant;
  end

  function Mdl = set.Constant(Mdl, value)
     validateattributes(value, {'double'}, {'scalar'}, '', 'Constant');
     Mdl.PrivateConstant = value;
     Mdl.FitInformation  = [];  % Indicate that the model is NOT estimated.
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function value = get.Variance(Mdl)
     value = Mdl.PrivateVariance;
  end

  function Mdl = set.Variance(Mdl, value)
     validateattributes(value, {'double' 'garch' 'gjr' 'egarch'}, {'scalar'}, '', 'Variance');
     updateDescription   = Mdl.Description == getModelSummary(Mdl);
     Mdl.PrivateVariance = value;
     try
       validateModel(Mdl, 'Variance', value);  % No model is returned (we only want error checking!)
     catch exception
       exception.throwAsCaller();
     end
     if ~isa(value, 'double')
        Mdl.PrivateDistribution = value.Distribution;
     end
     Mdl.FitInformation  = [];  % Indicate that the model is NOT estimated.
     if updateDescription       % Needed because conditional variance models update distribution too!
        Mdl.Description = getModelSummary(Mdl);
     end
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function value = get.Beta(Mdl)
     value = Mdl.PrivateBeta;
  end

  function Mdl = set.Beta(Mdl, value)
     if any(size(value))
        validateattributes(value, {'double'}, {'vector'}, '', 'Beta');
     else
        validateattributes(value, {'double'}, {}, '', 'Beta'); % Allow Mdl.Beta = [] for backward compatibility
     end
     updateDescription  = Mdl.Description == getModelSummary(Mdl);
     Mdl.PrivateBeta    = value(:)';    % Ensure a row vector
     Mdl.FitInformation = [];           % Indicate that the model is NOT estimated.
     if updateDescription
        Mdl.Description = getModelSummary(Mdl);
     end
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function value = get.D(Mdl)
     value = Mdl.PrivateD;
  end

  function Mdl = set.D(Mdl, value)
     validateattributes(value, {'double'}, {'scalar' 'nonnegative' 'integer'}, '', 'D');
     updateDescription  = Mdl.Description == getModelSummary(Mdl);
     Mdl.PrivateD       = value;
     Mdl.FitInformation = [];  % Indicate that the model is NOT estimated.
     if updateDescription
        Mdl.Description = getModelSummary(Mdl);
     end
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function value = get.Seasonality(Mdl)
     value = Mdl.PrivateSeasonality;
  end

  function Mdl = set.Seasonality(Mdl, value)
     validateattributes(value, {'double'}, {'scalar' 'nonnegative' 'integer'}, '', 'Seasonality')
     updateDescription      = Mdl.Description == getModelSummary(Mdl);
     Mdl.PrivateSeasonality = value;
     Mdl.FitInformation     = [];  % Indicate that the model is NOT estimated.
     if updateDescription
        Mdl.Description = getModelSummary(Mdl);
     end
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function coefficients = get.AR(Mdl)
     if Mdl.LagOpLHS{1}.Degree > 0
        coefficients = toCellArray(reflect(Mdl.LagOpLHS{1}));
        coefficients = coefficients(2:end);   % Coefficients at lags 1, 2, ...
     else
        coefficients = {};
     end
  end

  function Mdl = set.AR(Mdl, value)
     if isempty(value)
        validateattributes(value, {'cell'}, {}, '', 'AR');
     else
        validateattributes(value, {'cell'}, {'vector'}, '', 'AR');
     end
     updateDescription = Mdl.Description == getModelSummary(Mdl);
     for i = 1:numel(value)
         if isempty(value{i})
            value{i} = 0;
         else
            if ~isscalar(value{i})
               error(message('econ:arima:arima:NonScalarARCoefficients'))
            end
         end
     end
     Mdl.LagOpLHS{1} = LagOp([Mdl.LagOpLHS{1}.Coefficients{0} cellfun(@uminus, value(:)', 'uniformoutput', false)]);
     try
       validateModel(Mdl, 'AR', value);  % No model is returned (we only want error checking!)
     catch exception
       exception.throwAsCaller();
     end
     Mdl.FitInformation = [];  % Indicate that the model is NOT estimated.
     if updateDescription
        Mdl.Description = getModelSummary(Mdl);
     end
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function coefficients = get.SAR(Mdl)
     if Mdl.LagOpLHS{2}.Degree > 0
        coefficients = toCellArray(reflect(Mdl.LagOpLHS{2}));
        coefficients = coefficients(2:end);   % Coefficients at lags 1, 2, ...
     else
        coefficients = {};
     end
  end

  function Mdl = set.SAR(Mdl, value)
     if isempty(value)
        validateattributes(value, {'cell'}, {}, '', 'SAR');
     else
        validateattributes(value, {'cell'}, {'vector'}, '', 'SAR');
     end
     updateDescription = Mdl.Description == getModelSummary(Mdl);
     for i = 1:numel(value)
         if isempty(value{i})
            value{i} = 0;
         else
            if ~isscalar(value{i})
               error(message('econ:arima:arima:NonScalarSARCoefficients'))
            end
         end
     end
     Mdl.LagOpLHS{2} = LagOp([Mdl.LagOpLHS{2}.Coefficients{0} cellfun(@uminus, value(:)', 'uniformoutput', false)]);
     try
       validateModel(Mdl, 'SAR', value);  % No model is returned (we only want error checking!)
     catch exception
       exception.throwAsCaller();
     end
     Mdl.FitInformation = [];  % Indicate that the model is NOT estimated.
     if updateDescription
        Mdl.Description = getModelSummary(Mdl);
     end
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function coefficients = get.MA(Mdl)
     if Mdl.LagOpRHS{1}.Degree > 0
        coefficients = toCellArray(Mdl.LagOpRHS{1});
        coefficients = coefficients(2:end);   % Coefficients at lags 1, 2, ...
     else
        coefficients = {};
     end
  end

  function Mdl = set.MA(Mdl, value)
     if isempty(value)
        validateattributes(value, {'cell'}, {}, '', 'MA');
     else
        validateattributes(value, {'cell'}, {'vector'}, '', 'MA');
     end
     updateDescription = Mdl.Description == getModelSummary(Mdl);
     for i = 1:numel(value)
         if isempty(value{i})
            value{i} = 0;
         else
            if ~isscalar(value{i})
               error(message('econ:arima:arima:NonScalarMACoefficients'))
            end
         end
     end
     Mdl.LagOpRHS{1} = LagOp([Mdl.LagOpRHS{1}.Coefficients{0} value(:)']);
     try
       validateModel(Mdl, 'MA', value);  % No model is returned (we only want error checking!)
     catch exception
       exception.throwAsCaller();
     end
     Mdl.FitInformation = [];  % Indicate that the model is NOT estimated.
     if updateDescription
        Mdl.Description = getModelSummary(Mdl);
     end
  end
  
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function coefficients = get.SMA(Mdl)
     if Mdl.LagOpRHS{2}.Degree > 0
        coefficients = toCellArray(Mdl.LagOpRHS{2});
        coefficients = coefficients(2:end);   % Coefficients at lags 1, 2, ...
     else
        coefficients = {};
     end
  end

  function Mdl = set.SMA(Mdl, value)
     if isempty(value)
        validateattributes(value, {'cell'}, {}, '', 'SMA');
     else
        validateattributes(value, {'cell'}, {'vector'}, '', 'SMA');
     end
     updateDescription = Mdl.Description == getModelSummary(Mdl);
     for i = 1:numel(value)
         if isempty(value{i})
            value{i} = 0;
         else
            if ~isscalar(value{i})
               error(message('econ:arima:arima:NonScalarSMACoefficients'))
            end
         end
     end
     Mdl.LagOpRHS{2} = LagOp([Mdl.LagOpRHS{2}.Coefficients{0} value(:)']);
     try
       validateModel(Mdl, 'SMA', value);  % No model is returned (we only want error checking!)
     catch exception
       exception.throwAsCaller();
     end
     Mdl.FitInformation = [];  % Indicate that the model is NOT estimated.
     if updateDescription
        Mdl.Description = getModelSummary(Mdl);
     end
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function average = get.UnconditionalMean(Mdl)
     AR          = Mdl.getLagOp('Compound AR');
     isARstable  = isStable(AR);                 % Determine if the process is AR stable
     AR          = cell2mat(toCellArray(AR));    % Get AR coefficients as a vector
     nPredictors = size(Mdl.PrivateBeta,2);      % # of predictors
     if isARstable && (sum(AR) ~= 0) && (nPredictors == 0)
        try
           average = Mdl.PrivateConstant / sum(AR);
        catch
           average = nan;
        end
     else
        average = nan;
     end
  end

end  % METHODS Block


methods (Access = public)

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function Mdl = arima(varargin)
%ARIMA Construct an ARIMA(P,D,Q) model

if nargin == 0
% 
%  MATLAB classes should construct a scalar object in a default state in
%  the absence of input arguments. 
%
%  In this case, the default constructor syntax creates an ARIMA(0,0,0) model 
%  with a Gaussian conditional probability distribution, an undefined additive 
%  constant, and an undefined constant-variance model. 
%
%  Therefore, the default ARIMA model represents is time series of mean-zero
%  unit-variance residuals.
%
   Mdl.LagOpLHS    = {LagOp(1) LagOp(1)};
   Mdl.LagOpRHS    = {LagOp(1) LagOp(1)};
   Mdl.Constant    = NaN;
   Mdl.Variance    = NaN;
   Mdl.Description = getModelSummary(Mdl);
   return
end

%
% Validate input parameters.
%

if isnumeric(varargin{1})  && (nargin == 3)  % Is it the short-hand syntax?

%
%  Validate short-hand syntax, Mdl = ARIMA(P,D,Q).
%
   parser = inputParser;
   parser.addRequired('P', @(x) validateattributes(x, {'double'}, {'scalar' 'nonnegative' 'integer'}, '', 'P'));
   parser.addRequired('D', @(x) validateattributes(x, {'double'}, {'scalar' 'nonnegative' 'integer'}, '', 'D'));
   parser.addRequired('Q', @(x) validateattributes(x, {'double'}, {'scalar' 'nonnegative' 'integer'}, '', 'Q'));
   parser.parse(varargin{:});

   P = parser.Results.P;
   Q = parser.Results.Q;

   Mdl.LagOpLHS    = {LagOp([1 nan(1,P)]) LagOp(1)};
   Mdl.LagOpRHS    = {LagOp([1 nan(1,Q)]) LagOp(1)};
   Mdl.PrivateD    = parser.Results.D;
   Mdl.Description = getModelSummary(Mdl);

elseif ischar(varargin{1}) || ( isstring(varargin{1}) && isscalar(varargin{1}) )

%
%  Validate long-hand syntax of parameter name-value pairs.
%

   try
     Mdl = validateModel(Mdl, varargin{:});
   catch exception
     exception.throwAsCaller();
   end

else
   error(message('econ:arima:arima:InvalidInputSyntax'))
end

end % Constructor


end % Methods (Access = public)


methods (Access = protected)

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function displayScalarObject(Mdl)

%
% Override matlab.mixin.CustomDisplay.displayScalarObject method.
%

maxToPrint = 3;  % Maximum # of coefficients to explicitly print

disp(matlab.mixin.CustomDisplay.getSimpleHeader(Mdl));

spaces = '                    ';
fprintf([spaces(1:5) 'Description: "%s"\n'], Mdl.Description)

if strcmpi(Mdl.PrivateDistribution.Name, 'Gaussian')
   fprintf([spaces(1:4) 'Distribution: Name = "%s"\n'], Mdl.PrivateDistribution.Name)
else
   fprintf([spaces(1:4) 'Distribution: Name = "%s", DoF = %g\n'], Mdl.PrivateDistribution.Name, Mdl.PrivateDistribution.DoF)
end

fprintf([spaces(1:8) '       P: %d\n'], Mdl.P)
fprintf([spaces(1:8) '       D: %d\n'], Mdl.PrivateD)
fprintf([spaces(1:8) '       Q: %d\n'], Mdl.Q)

fprintf([spaces(1:8) 'Constant: %g\n'], Mdl.PrivateConstant)

% 
% Print AR coefficient information.
%

C = reflect(Mdl.LagOpLHS{1});
C = C.Coefficients;
L = Mdl.LagOpLHS{1}.Lags;
L = L(:,L > 0);
N = numel(L);                   % Get the number of included lags.

if N == 0
   fprintf([spaces(1:14) 'AR: {}\n'])
else
   format = repmat(' %g', 1, N);
   if N == 1
      fprintf([spaces(1:14) 'AR: {' format(2:end) '}' ' at lag [' format(2:end) ']\n'], C{L}, L)
   else
      fprintf([spaces(1:14) 'AR: {' format(2:end) '}' ' at lags [' format(2:end) ']\n'], C{L}, L)
   end
end

% 
% Print SAR coefficient information.
%

C = reflect(Mdl.LagOpLHS{2});
C = C.Coefficients;
L = Mdl.LagOpLHS{2}.Lags;
L = L(:,L > 0);
N = numel(L);                   % Get the number of included lags.

if N == 0
   fprintf([spaces(1:13) 'SAR: {}\n'])
else
   format = repmat(' %g', 1, N);
   if N == 1
      fprintf([spaces(1:13) 'SAR: {' format(2:end) '}' ' at lag [' format(2:end) ']\n'], C{L}, L)
   else
      fprintf([spaces(1:13) 'SAR: {' format(2:end) '}' ' at lags [' format(2:end) ']\n'], C{L}, L)
   end
end

% 
% Print MA coefficient information.
%

C = Mdl.LagOpRHS{1}.Coefficients;
L = Mdl.LagOpRHS{1}.Lags;
L = L(:,L > 0);
N = numel(L);                   % Get the number of included lags.

if N == 0
   fprintf([spaces(1:14) 'MA: {}\n'])
else
   format = repmat(' %g', 1, N);
   if N == 1
      fprintf([spaces(1:14) 'MA: {' format(2:end) '}' ' at lag [' format(2:end) ']\n'], C{L}, L)
   else
      fprintf([spaces(1:14) 'MA: {' format(2:end) '}' ' at lags [' format(2:end) ']\n'], C{L}, L)
   end
end

% 
% Print SMA coefficient information.
%

C = Mdl.LagOpRHS{2}.Coefficients;
L = Mdl.LagOpRHS{2}.Lags;
L = L(:,L > 0);
N = numel(L);                   % Get the number of included lags.

if N == 0
   fprintf([spaces(1:13) 'SMA: {}\n'])
else
   format = repmat(' %g', 1, N);
   if N == 1
      fprintf([spaces(1:13) 'SMA: {' format(2:end) '}' ' at lag [' format(2:end) ']\n'], C{L}, L)
   else
      fprintf([spaces(1:13) 'SMA: {' format(2:end) '}' ' at lags [' format(2:end) ']\n'], C{L}, L)
   end
end

%
% Print information about seasonality and regression coefficients.
%

fprintf([spaces(1:5) 'Seasonality: %d\n'], Mdl.PrivateSeasonality)

if numel(Mdl.PrivateBeta) == 0
   fprintf([spaces(1:12) 'Beta: [%d' char(215) '%d]\n'], size(Mdl.PrivateBeta))
else
   format = repmat(' %g', 1, numel(Mdl.PrivateBeta));
   fprintf([spaces(1:12) 'Beta: [' format(2:end) ']\n'], Mdl.PrivateBeta(:)')
end

if isa(Mdl.PrivateVariance, 'double')
   fprintf([spaces(1:8) 'Variance: %g\n'], Mdl.PrivateVariance)
else
   fprintf([spaces(1:8) 'Variance: [%s(%d,%d) Model]\n'], upper(class(Mdl.PrivateVariance)), Mdl.PrivateVariance.P, Mdl.PrivateVariance.Q)
end

%
% Print the model summary string if needed.
%

summary = getModelSummary(Mdl);

if Mdl.Description ~= summary
   disp(' ')
   fprintf("   " + '<strong>' + summary + '</strong>\n');
end

end % Display Scalar Object

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function Mdl = validateModel(Mdl, varargin)
%
% Validate an input model under the one of the following conditions:
%
%  o A variable-length list of parameter name-value pairs is passed when 
%    called from the constructor.
%
%  o A property SET method has been called to validate an assignment, but
%    no model is returned (we only want error checking!). There is only one
%    additional N-V pair passed (the assigned parameter of interest).
%
%  o Either the "vectorToModel" or "setLagOp" utility method has been called 
%    to validate a model. An output model is returned but no additional inputs 
%    are passed.
%

parser = inputParser;
parser.addParameter('Constant'    , Mdl.PrivateConstant    , @(x) validateattributes(x, {'double'}, {'scalar'}, '', 'Constant'));
parser.addParameter('D'           , Mdl.PrivateD           , @(x) validateattributes(x, {'double'}, {'scalar' 'nonnegative' 'integer'}, '', 'D'));
parser.addParameter('AR'          , Mdl.LagOpLHS{1}        , @(x) validateattributes(x, {'double' 'cell'}, {}, '', 'AR'));
parser.addParameter('SAR'         , Mdl.LagOpLHS{2}        , @(x) validateattributes(x, {'double' 'cell'}, {}, '', 'SAR'));
parser.addParameter('MA'          , Mdl.LagOpRHS{1}        , @(x) validateattributes(x, {'double' 'cell'}, {}, '', 'MA'));
parser.addParameter('SMA'         , Mdl.LagOpRHS{2}        , @(x) validateattributes(x, {'double' 'cell'}, {}, '', 'SMA'));
parser.addParameter('Distribution', Mdl.PrivateDistribution, @(x) validateattributes(x, {'char' 'string' 'struct'}, {}, '', 'Distribution'));
parser.addParameter('ARLags'      , []                     , @(x) validateattributes(x, {'double'}, {'integer'}, '', 'ARLags'));
parser.addParameter('MALags'      , []                     , @(x) validateattributes(x, {'double'}, {'integer'}, '', 'MALags'));
parser.addParameter('SARLags'     , []                     , @(x) validateattributes(x, {'double'}, {'integer'}, '', 'SARLags'));
parser.addParameter('SMALags'     , []                     , @(x) validateattributes(x, {'double'}, {'integer'}, '', 'SMALags'));
parser.addParameter('Seasonality' , Mdl.PrivateSeasonality , @(x) validateattributes(x, {'double'}, {'scalar' 'nonnegative' 'integer'}, '', 'Seasonality'));
parser.addParameter('Variance'    , Mdl.PrivateVariance    , @(x) validateattributes(x, {'double' 'garch' 'gjr' 'egarch'}, {'scalar'}, '', 'Variance'));
parser.addParameter('Beta'        , Mdl.PrivateBeta        , @(x) validateattributes(x, {'double'}, {}, '', 'Beta'));
parser.addParameter('Description' , ''                     , @(x) validateattributes(x, {'char' 'string'}         , {'scalartext'}       , '', 'Description'));

parser.parse(varargin{:});

AR           = parser.Results.AR;
SAR          = parser.Results.SAR;
MA           = parser.Results.MA;
SMA          = parser.Results.SMA;
ARLags       = parser.Results.ARLags;
MALags       = parser.Results.MALags;
SARLags      = parser.Results.SARLags;
SMALags      = parser.Results.SMALags;
variance     = parser.Results.Variance;
distribution = parser.Results.Distribution;
description  = string(parser.Results.Description);
  
if ~isempty(ARLags) && ~isvector(ARLags)
   error(message('econ:arima:arima:NonVectorARLags'))
end

if ~isempty(SARLags) && ~isvector(SARLags)
   error(message('econ:arima:arima:NonVectorSARLags'))       
end

if ~isempty(MALags) && ~isvector(MALags)
   error(message('econ:arima:arima:NonVectorMALags'))       
end

if ~isempty(SMALags) && ~isvector(SMALags)
   error(message('econ:arima:arima:NonVectorSMALags'))       
end

%
% Check the variance model.
%

if isa(variance, 'double')                                     % Check constant-variance models

   if variance <= 0
      error(message('econ:arima:arima:NonPositiveVariance'))
   end

elseif any(strcmp(class(variance), {'garch' 'gjr' 'egarch'}))  % Check conditional variance models

    if variance.Offset ~= 0
       error(message('econ:arima:arima:NonZeroVarianceOffset'))
    end
%
%   Enforce consistency between probability distributions.
%
%   The ARIMA model contains a conditional variance model, so set the
%   distribution of the variance model to that of the ARIMA model if explicitly 
%   specified by the user, allowing the outer, or container, ARIMA model to 
%   take precedence. 
%
%   However, if the user did NOT specify a distribution for the ARIMA model, 
%   then allow the distribution of the contained variance model to migrate to 
%   the ARIMA model.
%
    if ~any(strcmpi('Distribution', parser.UsingDefaults))     % Did the user specify a distribution?
       variance.Distribution = distribution;
    else
       distribution = variance.Distribution;
    end
    
else
   error(message('econ:arima:arima:InvalidVarianceType'))
end

%
% Check and assign model coefficients and corresponding lags.
%
% In the following, if the user specified coefficients without lags, then 
% derive lags from the coefficients. Similarly, if the user specified lags 
% without coefficients, then assign NaNs to coefficients consistent with 
% the number of lags.
%

%
% Additional checking for AR coefficients.
%

if isnumeric(AR)
   AR = num2cell(AR);
end

if isa(AR, 'LagOp')
   AR = toCellArray(reflect(AR));
   AR = AR(2:end);
end

if length(AR) == numel(AR)
   AR = AR(:).';
else
   error(message('econ:arima:arima:NonVectorARCoefficients'))
end

for i = 1:numel(AR)
     if ~isscalar(AR{i})
        error(message('econ:arima:arima:NonScalarARCoefficients'))
     end
end

if any(strcmpi('ARLags', parser.UsingDefaults))
   ARLags = 0:numel(AR);
else
   if any(strcmpi('AR', parser.UsingDefaults))
      AR = num2cell(nan(1,numel(ARLags)));
   end

   if numel(ARLags) ~= numel(AR)
      error(message('econ:arima:arima:InconsistentARSpec'))
   end

   Lags = unique(ARLags, 'first');

   if any(ARLags == 0) || (numel(ARLags) ~= numel(Lags))
      error(message('econ:arima:arima:InvalidARLags'))
   end
   ARLags = [0 ARLags];
end

%
% Additional checking for SAR coefficients.
%

if isnumeric(SAR)
   SAR = num2cell(SAR);
end

if isa(SAR, 'LagOp')
   SAR = toCellArray(reflect(SAR));
   SAR = SAR(2:end);
end

if length(SAR) == numel(SAR)
   SAR = SAR(:).';
else
   error(message('econ:arima:arima:NonVectorSARCoefficients'))
end

for i = 1:numel(SAR)
    if ~isscalar(SAR{i})
       error(message('econ:arima:arima:NonScalarSARCoefficients'))
    end
end

if any(strcmpi('SARLags', parser.UsingDefaults))
   SARLags = 0:numel(SAR);
else
   if any(strcmpi('SAR', parser.UsingDefaults))
      SAR = num2cell(nan(1,numel(SARLags)));
   end

   if numel(SARLags) ~= numel(SAR)
      error(message('econ:arima:arima:InconsistentSARSpec'))
   end

   Lags = unique(SARLags, 'first');

   if any(SARLags == 0) || (numel(SARLags) ~= numel(Lags))
      error(message('econ:arima:arima:InvalidSARLags'))
   end
   SARLags = [0 SARLags];
end

%
% Additional checking for MA coefficients.
%

if isnumeric(MA)
   MA = num2cell(MA);
end

if isa(MA, 'LagOp')
   MA = toCellArray(MA);
   MA = MA(2:end);
end

if length(MA) == numel(MA)
   MA = MA(:).';
else
   error(message('econ:arima:arima:NonVectorMACoefficients'))
end

for i = 1:numel(MA)
    if ~isscalar(MA{i})
       error(message('econ:arima:arima:NonScalarMACoefficients'))
    end
end   

if any(strcmpi('MALags', parser.UsingDefaults))
   MALags = 0:numel(MA);
else
   if any(strcmpi('MA', parser.UsingDefaults))
      MA = num2cell(nan(1,numel(MALags)));
   end

   if numel(MALags) ~= numel(MA)
      error(message('econ:arima:arima:InconsistentMASpec'))
   end

   Lags = unique(MALags);

   if any(MALags == 0) || (numel(MALags) ~= numel(Lags))
      error(message('econ:arima:arima:InvalidMALags'))
   end
   MALags = [0 MALags];
end

%
% Additional checking for SMA coefficients.
%

if isnumeric(SMA)
   SMA = num2cell(SMA);
end

if isa(SMA, 'LagOp')
   SMA = toCellArray(SMA);
   SMA = SMA(2:end);
end

if length(SMA) == numel(SMA)
   SMA = SMA(:).';
else
   error(message('econ:arima:arima:NonVectorSMACoefficients'))
end

for i = 1:numel(SMA)
    if ~isscalar(SMA{i})
       error(message('econ:arima:arima:NonScalarSMACoefficients'))
    end
end   

if any(strcmpi('SMALags', parser.UsingDefaults))
   SMALags = 0:numel(SMA);
else
   if any(strcmpi('SMA', parser.UsingDefaults))
      SMA = num2cell(nan(1,numel(SMALags)));
   end

   if numel(SMALags) ~= numel(SMA)
      error(message('econ:arima:arima:InconsistentSMASpec'))
   end

   Lags = unique(SMALags);

   if any(SMALags == 0) || (numel(SMALags) ~= numel(Lags))
      error(message('econ:arima:arima:InvalidSMALags'))
   end
   SMALags = [0 SMALags];
end

%
% No output model is an indication that this validation method is called ONLY 
% to perform error checking, and so updating the model is unnecessary.
%
% If an updated model is returned, then create the underlying Lag Operator 
% Polynomials and assign or auto-generate the model description.
%
% Additional properties are assigned using PRIVATE interface to avoid redundant 
% error checking.

if nargout > 0
%
%  Create the underlying Lag Operator Polynomials.
%
   Mdl.LagOpLHS{1} = LagOp([1 cellfun(@uminus, AR , 'uniformoutput', false)], 'Lags', ARLags);
   Mdl.LagOpLHS{2} = LagOp([1 cellfun(@uminus, SAR, 'uniformoutput', false)], 'Lags', SARLags);
   Mdl.LagOpRHS{1} = LagOp([1  MA], 'Lags',  MALags);
   Mdl.LagOpRHS{2} = LagOp([1 SMA], 'Lags', SMALags);
%   
%  Assign using PRIVATE interface to avoid redundant error checking.
%
   Mdl.PrivateConstant    = parser.Results.Constant;
   Mdl.PrivateD           = parser.Results.D;
   Mdl.PrivateSeasonality = parser.Results.Seasonality;
   Mdl.PrivateVariance    = variance;
%
%  Assign or auto-generate the model description.
%
   if ~any(strcmpi('Description', parser.UsingDefaults))
      Mdl.Description = description;
   else
      Mdl.Description = getModelSummary(Mdl);
   end
end

%
% Assign remaining properties using the PUBLIC interface to perform
% property-specific error checking NOT performed in this validation method.
%

Mdl.Beta         = parser.Results.Beta;
Mdl.Distribution = distribution;

%
%  Indicate stationarity/invertibility constraints for coefficients.
%

if all(~isnan(cell2mat(AR))) && ~isempty(cell2mat(AR)) && ~isStable(Mdl.LagOpLHS{1})
   error(message('econ:arima:arima:UnstableAR'))
end

if all(~isnan(cell2mat(SAR))) && ~isempty(cell2mat(SAR)) && ~isStable(Mdl.LagOpLHS{2})
   error(message('econ:arima:arima:UnstableSAR'))
end

if all(~isnan(cell2mat(MA))) && ~isempty(cell2mat(MA)) && ~isStable(Mdl.LagOpRHS{1})
   error(message('econ:arima:arima:UnstableMA'))
end

if all(~isnan(cell2mat(SMA))) && ~isempty(cell2mat(SMA)) && ~isStable(Mdl.LagOpRHS{2})
   error(message('econ:arima:arima:UnstableSMA'))
end


end % Validate Model


end % Methods (Protected)


methods (Hidden)

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function varargout = getModelSummary(Mdl)
%
% Create the model summary string used by the "displayScalarObject" 
% and "summarize" methods.
%

summary = "   ";

if isempty(Mdl.PrivateBeta) || all(Mdl.PrivateBeta(:) == 0)
   summary = summary + "ARIMA";
else
   summary = summary + "ARIMAX";
end

if Mdl.PrivateSeasonality > 0
   summary = summary + "(%d,%d,%d) Model Seasonally Integrated";
else
   summary = summary + "(%d,%d,%d) Model";
end

if (Mdl.LagOpLHS{2}.Degree > 0) && (Mdl.LagOpRHS{2}.Degree > 0)
   summary = summary + " with Seasonal AR(%d) and MA(%d)";
   summary = sprintf(summary, Mdl.LagOpLHS{1}.Degree, Mdl.PrivateD, Mdl.LagOpRHS{1}.Degree, Mdl.LagOpLHS{2}.Degree, Mdl.LagOpRHS{2}.Degree);
elseif Mdl.LagOpLHS{2}.Degree > 0
   summary = summary + " with Seasonal AR(%d)";
   summary = sprintf(summary, Mdl.LagOpLHS{1}.Degree, Mdl.PrivateD, Mdl.LagOpRHS{1}.Degree, Mdl.LagOpLHS{2}.Degree);
elseif Mdl.LagOpRHS{2}.Degree > 0
   summary = summary + " with Seasonal MA(%d)";
   summary = sprintf(summary, Mdl.LagOpLHS{1}.Degree, Mdl.PrivateD, Mdl.LagOpRHS{1}.Degree, Mdl.LagOpRHS{2}.Degree);
else
   summary = sprintf(summary, Mdl.LagOpLHS{1}.Degree, Mdl.PrivateD, Mdl.LagOpRHS{1}.Degree);
end

if strcmpi(Mdl.PrivateDistribution.Name, 'T')    % Is it a t distribution?
   summary = summary + " (t Distribution)";
else
   summary = summary + " (Gaussian Distribution)";
end

%
% Create & print the model summary string.
%

if nargout > 0
   varargout = {strtrim(summary)};
else
   fprintf("<strong>" + summary + "</strong>\n");
end

end % Get Model Summary

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function Mdl = vectorToModel(Mdl, coefficients, covariance)
%VECTORTOMODEL Pack estimated parameters into an ARIMA model
%
% Syntax:
%
%   Mdl = vectorToModel(Mdl, coefficients, covariance)
%
% Description:
%
%   Given a vector of estimated coefficients, as output from FMINCON, pack
%   the coefficients into an ARIMA model. Also, compute the standard
%   errors from the estimated parameter error covariance matrix and store
%   the standard errors and related "fit information" in the output model.
%
% Input Arguments:
%
%   Mdl - ARIMA model specification to update with parameter estimates.
%
%   coefficients - A vector of estimated coefficients to pack into the model Mdl.
%
%   covariance - Variance-covariance matrix associated with model parameters 
%     known to the optimizer. The rows and columns associated with any 
%     parameters estimated by maximum likelihood contain the covariances of 
%     the estimation errors; the standard errors of the parameter estimates 
%     are the square root of the entries along the main diagonal. The rows 
%     and columns associated with any parameters held fixed as equality 
%     constraints contain zeros. 
%
% Output Arguments:
%
%   Mdl - An updated ARIMA model specification containing the parameter 
%     estimates and related fit information.
%

if iscell(coefficients)
   coefficients = cell2mat(coefficients);
end

coefficients = coefficients(:);      % Guarantee a column vector

updateDescription    =  Mdl.Description == getModelSummary(Mdl);                % Should the description be updated?
isVarianceConstant   = ~any(strcmp(class(Mdl.PrivateVariance), {'garch' 'gjr' 'egarch'})); % Is it a constant variance model?
isDistributionT      =  strcmpi(Mdl.PrivateDistribution.Name, 'T');             % Is it a t distribution?
isRegressionIncluded = ~isempty(Mdl.PrivateBeta);

%
% Determine the non-zero lags and corresponding coefficient indices.
%

LagsAR  = Mdl.LagOpLHS{1}.Lags;       % Lags of non-zero  AR coefficients
LagsMA  = Mdl.LagOpRHS{1}.Lags;       % Lags of non-zero  MA coefficients
LagsSAR = Mdl.LagOpLHS{2}.Lags;       % Lags of non-zero SAR coefficients
LagsSMA = Mdl.LagOpRHS{2}.Lags;       % Lags of non-zero SMA coefficients

LagsAR  = LagsAR(LagsAR > 0);         % Retain only positive lags
LagsMA  = LagsMA(LagsMA > 0);
LagsSAR = LagsSAR(LagsSAR > 0);
LagsSMA = LagsSMA(LagsSMA > 0);

nCoefficients = 1 + numel(LagsAR) + numel(LagsSAR) + ...
                    numel(LagsMA) + numel(LagsSMA) + size(Mdl.PrivateBeta,2);

iAR   = 2:(2 + numel(LagsAR) - 1);
iSAR  = (numel(iAR) + 2):(numel(iAR) + numel(LagsSAR) + 1);
iMA   = (numel(iAR) + numel(iSAR) + 2):(numel(iAR) + numel(iSAR) + numel(LagsMA) + 1);
iSMA  = (numel(iAR) + numel(iSAR) + numel(iMA) + 2):(numel(iAR) + numel(iSAR) + numel(LagsMA) + numel(LagsSMA) + 1);
iBeta = (numel(iAR) + numel(iSAR) + numel(iMA) + numel(iSMA) + 2):nCoefficients;

%
% Determine which parameters are estimated (TRUE) and which are not (FALSE).
%

isEstimated.Constant = isnan(Mdl.PrivateConstant);

if Mdl.LagOpLHS{1}.Degree > 0
   for i = 1:Mdl.LagOpLHS{1}.Degree
       isEstimated.AR{i} = isnan(Mdl.AR{i});
   end
else
   isEstimated.AR = {};
end

if Mdl.LagOpLHS{2}.Degree > 0
   for i = 1:Mdl.LagOpLHS{2}.Degree
       isEstimated.SAR{i} = isnan(Mdl.SAR{i});
   end
else
   isEstimated.SAR = {};
end

if Mdl.LagOpRHS{1}.Degree > 0
   for i = 1:Mdl.LagOpRHS{1}.Degree
       isEstimated.MA{i} = isnan(Mdl.MA{i});
   end
else
   isEstimated.MA = {};
end

if Mdl.LagOpRHS{2}.Degree > 0
   for i = 1:Mdl.LagOpRHS{2}.Degree
       isEstimated.SMA{i} = isnan(Mdl.SMA{i});
   end
else
   isEstimated.SMA = {};
end

%
% Unpack the coefficient vector and assign corresponding standard errors and
% update the reported ARIMA object.
%

sigma                   = sqrt(diag(covariance));  % Extract standard errors
Mdl.PrivateConstant     = coefficients(1);
standardErrors.Constant = sigma(1);

if isRegressionIncluded  % Is there a regression component?
   isEstimated.Beta    = isnan(Mdl.PrivateBeta);
   beta                = coefficients(iBeta);      % Regression coefficients
   Mdl.PrivateBeta     = beta(:)';                 % Ensure a row vector 
   standardErrors.Beta = sigma(iBeta)';
end

if isVarianceConstant
   isEstimated.Variance    = isnan(Mdl.PrivateVariance);
   Mdl.PrivateVariance     = coefficients(end - isDistributionT);
   standardErrors.Variance = sigma(end - isDistributionT);
else
   Mdl.PrivateVariance = vectorToModel(Mdl.PrivateVariance, coefficients((nCoefficients + 1):end), ...
                                                            covariance((nCoefficients + 1):end,(nCoefficients + 1):end));
end

if isDistributionT       % Is it a t distribution?
   isEstimated.DoF             = isnan(Mdl.PrivateDistribution.DoF);
   Mdl.PrivateDistribution.DoF = coefficients(end);  % DoF
   standardErrors.DoF          = sigma(end);
end

%
% Assign coefficients and their standard errors. 
%

Mdl.LagOpLHS{1}  = LagOp([1 -coefficients(iAR)' ], 'Lags', [0 LagsAR]);
Mdl.LagOpRHS{1}  = LagOp([0  coefficients(iMA)' ], 'Lags', [0 LagsMA]);
Mdl.LagOpLHS{2}  = LagOp([1 -coefficients(iSAR)'], 'Lags', [0 LagsSAR]);
Mdl.LagOpRHS{2}  = LagOp([0  coefficients(iSMA)'], 'Lags', [0 LagsSMA]);

if Mdl.LagOpLHS{1}.Degree > 0
   for i = 1:Mdl.LagOpLHS{1}.Degree
       if any(i == LagsAR)
          iLag = (i == LagsAR);
          standardErrors.AR{i} = sigma(iAR(iLag));
       else
          standardErrors.AR{i} = 0;
       end
   end
else
   standardErrors.AR = {};
end

if Mdl.LagOpLHS{2}.Degree > 0
   for i = 1:Mdl.LagOpLHS{2}.Degree
       if any(i == LagsSAR)
          iLag = (i == LagsSAR);
          standardErrors.SAR{i} = sigma(iSAR(iLag));
       else
          standardErrors.SAR{i} = 0;
       end
   end
else
   standardErrors.SAR = {};
end

if Mdl.LagOpRHS{1}.Degree > 0
   for i = 1:Mdl.LagOpRHS{1}.Degree
       if any(i == LagsMA)
          iLag = (i == LagsMA);
          standardErrors.MA{i} = sigma(iMA(iLag));
       else
          standardErrors.MA{i} = 0;
       end
   end
else
   standardErrors.MA = {};
end

if Mdl.LagOpRHS{2}.Degree > 0
   for i = 1:Mdl.LagOpRHS{2}.Degree
       if any(i == LagsSMA)
          iLag = (i == LagsSMA);
          standardErrors.SMA{i} = sigma(iSMA(iLag));
       else
          standardErrors.SMA{i} = 0;
       end
   end
else
   standardErrors.SMA = {};
end

%
% Validate the updated model and assign post-estimation information.
%

variance            = Mdl.PrivateVariance;  % Save the private FitInformation
Mdl                 = validateModel(Mdl);
Mdl.PrivateVariance = variance;             % Restore the private FitInformation

Mdl.FitInformation.IsEstimated    = isEstimated;
Mdl.FitInformation.StandardErrors = standardErrors;

if updateDescription              % Update model description if NOT user-specified
   Mdl.Description = getModelSummary(Mdl);
end

end % Vector-to-model

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function polynomial = getLagOp(Mdl, name)
%GETLAGOP Get lag operator polynomials of ARIMA models.
%
% Syntax:
%
%   polynomial = getLagOp(Mdl,name)
%
% Description:
%
%   Given the name of a component polynomial of the input model, return the
%   underlying lag operator polynomial (LagOp object).
%
% Input Arguments:
%
%   Mdl - The model whose component polynomial is requested.
%
%   name - The name of the polynomial (case-sensitive character string).
%     The available options are 'AR' (non-seasonal autoregressive), 'SAR' 
%     (seasonal autoregressive), 'MA' (non-seasonal moving average), 'SMA' 
%     (seasonal moving average), 'Integrated Non-Seasonal' (polynomial
%     associated with D, the degree of non-seasonal differencing),
%     'Integrated Seasonal' (the polynomial associated with Seasonality, 
%     the degree of non-seasonal differencing), 'Compound AR' (the compound
%     autoregressive polynomial, including all seasonal and non-seasonal 
%     components), and 'Compound MA' (the compound moving average polynomial, 
%     including all seasonal and non-seasonal components).
%
% Output Arguments:
%
%   polynomial - A lag operator polynomial associated with the input name.

   switch name                           % Name of the component polynomial

      case {'AR'}                        % AR  (Non-Seasonal Autoregressive)
         polynomial = Mdl.LagOpLHS{1};

      case {'SAR'}                       % SAR (Seasonal Autoregressive)
         polynomial = Mdl.LagOpLHS{2};

      case {'MA'}                        % MA  (Non-Seasonal Moving Average)
         polynomial = Mdl.LagOpRHS{1};

      case {'SMA'}                       % SMA (Seasonal Moving Average)
         polynomial = Mdl.LagOpRHS{2};

      case 'Integrated Non-Seasonal'
         D = Mdl.PrivateD;
         if D > 0
            polynomial = LagOp([1 -1]);
            for i = 2:D
                polynomial = polynomial * LagOp([1 -1]);
            end
         else
            polynomial = LagOp(1);
         end

      case 'Integrated Seasonal'
         S = Mdl.PrivateSeasonality;
         if S > 0
            polynomial = LagOp([1 -1], 'Lags', [0 S]);
         else
            polynomial = LagOp(1);
         end

      case 'Compound AR'
         polynomial = Mdl.LagOpLHS{1} * Mdl.getLagOp('Integrated Seasonal') * ...
                 Mdl.LagOpLHS{2} * Mdl.getLagOp('Integrated Non-Seasonal');

      case 'Compound MA'
         polynomial = Mdl.LagOpRHS{1} * Mdl.LagOpRHS{2};

      otherwise
         error(message('econ:arima:arima:InvalidPolynomialReference'))
   end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function Mdl = setLagOp(Mdl, name, polynomial, ~)
%SETLAGOP Set lag operator polynomials of ARIMA models.
%
% Syntax:
%
%   Mdl = setLagOp(Mdl,name,polynomial)
%
% Description:
%
%   Given the name of a component polynomial of the input model, update the
%   underlying lag operator polynomial (LagOp object).
%
% Input Arguments:
%
%   Mdl - The model whose component polynomial is updated.
%
%   name - The name of the polynomial (case-sensitive character string).
%     The available options are: 'AR' (non-seasonal autoregressive), 'SAR' 
%     (seasonal autoregressive), 'MA' (non-seasonal moving average), and 
%     'SMA' (seasonal moving average).
%
%   polynomial - A lag operator polynomial associated with the input name.
%
% Output Arguments:
%
%   Mdl - The model whose component polynomial is updated.

   switch name                               % Name of the component polynomial

      case {'AR'}                            % AR polynomial       
         Mdl.LagOpLHS{1} = polynomial;

      case {'SAR'}                           % SAR polynomial
         Mdl.LagOpLHS{2} = polynomial;

      case {'MA'}                            % MA polynomial
         Mdl.LagOpRHS{1} = polynomial;

      case {'SMA'}                           % SMA polynomial
         Mdl.LagOpRHS{2} = polynomial;

      otherwise
         error(message('econ:arima:arima:InvalidPolynomialReference'))
   end

   if nargin < 4
      Mdl = validateModel(Mdl);
   end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function print(Mdl, covariance)
%PRINT Print ARIMA model estimation results
%
% Syntax:
%
%   print(Mdl, VarCov)
%
% Description:
%
%   Given an ARIMA model containing parameter estimates and corresponding 
%   estimation error variance-covariance matrix, print a table of parameter
%   estimates, standard errors, and t-statistics.
%
% Input Arguments:
%
%   Mdl - ARIMA model specification.
%
%   VarCov - Estimation error variance-covariance matrix. VarCov is a square
%     matrix with a row and column associated with each parameter known to 
%     the optimizer. Known parameters include all parameters estimated as 
%     well as all coefficients held fixed throughout the optimization. 
%     Coefficients of lag operator polynomials at lags excluded from the 
%     model (held fixed at zero) are omitted from VarCov.
%
% Note:
%
%   o The parameters known to FMINCON and included in VarCov are ordered as 
%     follows:
%
%       - Constant
%       - Non-zero  AR coefficients at positive lags
%       - Non-zero SAR coefficients at positive lags
%       - Non-zero  MA coefficients at positive lags
%       - Non-zero SMA coefficients at positive lags
%       - Regression coefficients (models with regression components only)
%       - Variance parameters (scalar for constant-variance models, vector
%         of additional parameters otherwise)
%       - Degrees-of-freedom (t distributions only)
%

warning(message('econ:arima:arima:UseSummarizeInstead'))

if size(covariance,1) ~= size(covariance,2)
   error(message('econ:arima:arima:NonSquareCovarianceMatrix'))
end

%
% Determine which parameters are estimated and which are held fixed. 
%
% Parameters held fixed throughout estimation are associated with all-zero
% rows and columns in the error covariance matrix.
%

solve = sum(covariance,2) ~= 0;     % TRUE = estimated, FALSE = fixed

%
% Get model parameters and ensure coefficients are specified.
%

constant = Mdl.Constant;            % Model constant
beta     = Mdl.Beta;                % Regression coefficients
beta     = beta(:)';                % Guarantee a column vector

AR  = reflect(getLagOp(Mdl, 'AR'));
SAR = reflect(getLagOp(Mdl, 'SAR'));
MA  = getLagOp(Mdl, 'MA');
SMA = getLagOp(Mdl, 'SMA');

LagsAR  = AR.Lags;                  % Non-zero lags of non-seasonal AR coefficients
LagsSAR = SAR.Lags;                 % Non-zero lags of seasonal AR coefficients
LagsMA  = MA.Lags;                  % Non-zero lags of non-seasonal MA coefficients
LagsSMA = SMA.Lags;                 % Non-zero lags of seasonal MA coefficients

LagsAR  = LagsAR(LagsAR > 0);       % Retain only positive lags
LagsSAR = LagsSAR(LagsSAR > 0);
LagsMA  = LagsMA(LagsMA > 0);
LagsSMA = LagsSMA(LagsSMA > 0);

if isempty(LagsAR)
   AR = [];
else
   AR = AR.Coefficients;
   AR = [AR{LagsAR}];               % Non-zero AR coefficients (vector)
end

if isempty(LagsSAR)
   SAR = []; 
else
   SAR = SAR.Coefficients;
   SAR = [SAR{LagsSAR}];            % Non-zero SAR coefficients (vector)
end

if isempty(LagsMA)
   MA = [];
else
   MA = MA.Coefficients;
   MA = [MA{LagsMA}];               % Non-zero MA coefficients (vector)
end

if isempty(LagsSMA)
   SMA = []; 
else
   SMA = SMA.Coefficients;
   SMA = [SMA{LagsSMA}];            % Non-zero SMA coefficients (vector)
end

isDistributionT    =  strcmpi(Mdl.Distribution.Name, 'T');
isVarianceConstant = ~any(strcmp(class(Mdl.PrivateVariance), {'garch' 'gjr' 'egarch'}));  % Is it a constant variance model?

%
% Determine the number of coefficients associated with the ARIMA model (i.e.,
% excluding any parameters found in Mdl associated with the variance model
% and the degrees-of-freedom of t distributions).
%

nARIMA = 1 + numel(LagsAR) + numel(LagsSAR) + ...
             numel(LagsMA) + numel(LagsSMA) + numel(beta);

%
% Pack the vector of parameters associated with the ARIMA model; each element 
% of the parameters vector will be printed by THIS method. 
%
% In the event a t distribution is encountered, notice that the last element 
% of the parameters vector is always the degrees-of-freedom regardless of 
% whether the variance model is constant or conditional.
%

if isVarianceConstant
%
%  For constant-variance models, the number of elements in PARAMETERS will
%  equal nARIMA + 1 for Gaussian distributions and nARIMA + 2 for t distributions.
%   
   if isDistributionT
      parameters          = zeros(nARIMA + 2,1);
      parameters(end - 1) = Mdl.Variance;
      parameters(end)     = Mdl.Distribution.DoF;
   else
      parameters          = zeros(nARIMA + 1,1);
      parameters(end)     = Mdl.Variance;
   end
   
   parameters(1:nARIMA) = [constant AR  SAR  MA  SMA  beta];

   if numel(parameters) ~= numel(solve)
      error(message('econ:arima:arima:ModelCovarianceInconsistency'))
   end

else  % Conditional variance model
%
%  For conditional-variance models, the number of elements in PARAMETERS 
%  will be fewer than the number of elements in the SOLVE indicator because
%  the coefficients associated with the variance model are printed by the
%  print method of the variance model, and are therefore excluded from 
%  PARAMETERS.
%

   if isDistributionT
      parameters(nARIMA + 1) = Mdl.Distribution.DoF;
      parameters(1:nARIMA)   = [constant AR  SAR  MA  SMA  beta];
   else
      parameters = [constant AR  SAR  MA  SMA  beta];
   end

end

%
% Display summary information.
%

disp(' ')

if isempty(Mdl.Beta)
   name = 'ARIMA';
else
   name = 'ARIMAX';
end

if Mdl.Seasonality > 0
   s = ['    ' name '(%d,%d,%d) Model Seasonally Integrated'];
else
   s = ['    ' name '(%d,%d,%d) Model'];
end

if (Mdl.LagOpLHS{2}.Degree > 0) && (Mdl.LagOpRHS{2}.Degree > 0)
   s = [s ' with Seasonal AR(%d) and MA(%d):\n'];
   fprintf(s, Mdl.LagOpLHS{1}.Degree, Mdl.D, Mdl.LagOpRHS{1}.Degree, ...
              Mdl.LagOpLHS{2}.Degree, Mdl.LagOpRHS{2}.Degree)
   fprintf('    %s\n',  repmat('-', 1, numel(s) - 8))
elseif Mdl.LagOpLHS{2}.Degree > 0
   s = [s ' with Seasonal AR(%d):\n'];
   fprintf(s, Mdl.LagOpLHS{1}.Degree, Mdl.D, Mdl.LagOpRHS{1}.Degree, ...
              Mdl.LagOpLHS{2}.Degree)
   fprintf('    %s\n',  repmat('-', 1, numel(s) - 8))
elseif Mdl.LagOpRHS{2}.Degree > 0
   s = [s ' with Seasonal MA(%d):\n'];
   fprintf(s, Mdl.LagOpLHS{1}.Degree, Mdl.D, Mdl.LagOpRHS{1}.Degree, ...
              Mdl.LagOpRHS{2}.Degree)
   fprintf('    %s\n',  repmat('-', 1, numel(s) - 8))
else
   s = [s ':\n'];
   fprintf(s, Mdl.LagOpLHS{1}.Degree, Mdl.D, Mdl.LagOpRHS{1}.Degree)
   fprintf('    %s\n',  repmat('-', 1, numel(s) - 8))
end

fprintf('    Conditional Probability Distribution: %s\n\n'   , Mdl.Distribution.Name);

header =  ['                                  Standard          t     ' ;
           '     Parameter       Value          Error       Statistic ' ;
           '    -----------   -----------   ------------   -----------'];

disp(header)

%
% Format annotation strings and display the information. Standard errors
% and t-statistics held fixed by equality constraints during estimation
% have zero estimation error, and are designated 'Fixed'.
%

Fix    = ~solve;
errors = sqrt(diag(covariance));

if Fix(1)
   fprintf('     Constant   %12.6g  %12s   %12s\n', parameters(1), 'Fixed', 'Fixed')
else
   tStatistic = parameters(1) / errors(1); % Compute t-statistic for this parameter
   fprintf('     Constant   %12.6g  %12.6g   %12.6g\n', parameters(1), errors(1), tStatistic)
end

row = 2;

for i = LagsAR
    s = ['       ', repmat(' ', 1, (0 <= i) && (i < 10))];
    if Fix(row)
      fprintf('%sAR{%d}   %12.6g  %12s   %12s \n', s, i, parameters(row), 'Fixed', 'Fixed')
    else
       tStatistic = parameters(row) / errors(row);
       fprintf('%sAR{%d}   %12.6g  %12.6g   %12.6g\n', s, i, parameters(row), errors(row), tStatistic)
    end
    row = row + 1;
end

for i = LagsSAR
    s = ['      ', repmat(' ', 1, (0 <= i) && (i < 10))];
    if Fix(row)
       fprintf('%sSAR{%d}   %12.6g  %12s   %12s\n', s, i, parameters(row), 'Fixed', 'Fixed')
    else
       tStatistic = parameters(row) / errors(row);
       fprintf('%sSAR{%d}   %12.6g  %12.6g   %12.6g\n', s, i, parameters(row), errors(row), tStatistic)
    end
    row = row + 1;
end

for i = LagsMA
    s = ['       ', repmat(' ', 1, (0 <= i) && (i < 10))];
    if Fix(row)
       fprintf('%sMA{%d}   %12.6g  %12s   %12s \n', s, i, parameters(row), 'Fixed', 'Fixed')
    else
       tStatistic = parameters(row) / errors(row);
       fprintf('%sMA{%d}   %12.6g  %12.6g   %12.6g\n', s, i, parameters(row), errors(row), tStatistic)
    end
    row = row + 1;
end

for i = LagsSMA
    s = ['      ', repmat(' ', 1, (0 <= i) && (i < 10))];
    if Fix(row)
       fprintf('%sSMA{%d}   %12.6g  %12s   %12s\n', s, i, parameters(row), 'Fixed', 'Fixed')
    else
       tStatistic = parameters(row) / errors(row);
       fprintf('%sSMA{%d}   %12.6g  %12.6g   %12.6g\n', s, i, parameters(row), errors(row), tStatistic)
    end
    row = row + 1;
end

for i = 1:numel(beta)
    s = ['     ', repmat(' ', 1, (0 <= i) && (i < 10))];
    if Fix(row)
       fprintf('%sBeta(%d)   %12.6g  %12s   %12s\n', s, i, parameters(row), 'Fixed', 'Fixed')
    else
       tStatistic = parameters(row) / errors(row);
       fprintf('%sBeta(%d)   %12.6g  %12.6g   %12.6g\n', s, i, parameters(row), errors(row), tStatistic)
    end
    row = row + 1;
end

if isVarianceConstant
 
   if Fix(row)
      fprintf('     Variance   %12.6g  %12s   %12s\n', parameters(row), 'Fixed', 'Fixed')
   else
      tStatistic = parameters(row) / errors(row);
      fprintf('     Variance   %12.6g  %12.6g   %12.6g\n', parameters(row), errors(row), tStatistic)
   end

   if isDistributionT
      if Fix(end)
         fprintf('          DoF   %12.6g  %12s   %12s\n', parameters(end), 'Fixed', 'Fixed')
      else
         tStatistic = parameters(end) / errors(end);
         fprintf('          DoF   %12.6g  %12.6g   %12.6g\n', parameters(end), errors(end), tStatistic)
      end
   end

else  % Conditional variance model
 
   if isDistributionT
      if Fix(end)
         fprintf('          DoF   %12.6g  %12s   %12s\n', parameters(end), 'Fixed', 'Fixed')
      else
         tStatistic = parameters(end) / errors(end);
         fprintf('          DoF   %12.6g  %12.6g   %12.6g\n', parameters(end), errors(end), tStatistic)
      end
   end
 
   disp(' ')
   try
     Mdl.Variance.print(covariance((nARIMA + 1):end,(nARIMA + 1):end))
   catch exception
     exception.throwAsCaller();
   end
end

end % Print Method


end % Methods (Hidden)


methods (Static, Hidden) % Utilities

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function [nLogL, nLogLikelihoods, E] = nLogLikeGaussian(coefficients, Y, X, E, V, LagsAR, LagsMA, numPaths, maxPQ, T)

%
% Infer the residuals.
%

if isempty(X)

   E = internal.econ.arimaxMex(LagsAR, LagsMA, coefficients, maxPQ, Y', E')';
%
%  Note:
%
%  To improve runtime performance, the following code segment has been replaced 
%  by a private/undocumented C-MEX file called above, which does exactly the 
%  same thing. To restore the original MATLAB code, simply uncomment the 
%  following code segment:


%    I = ones(numPaths,1);
%
%    for t = (maxPQ + 1):T
%        data   = [I  Y(:,t - LagsAR)  E(:,t - LagsMA)];
%        E(:,t) = data * coefficients;
%    end

else
 
   arimaCoefficients   = coefficients(1:(end - size(X,1)));
   regressCoefficients = coefficients((end - size(X,1) + 1):end)';

   E = internal.econ.arimaxMex(LagsAR, LagsMA, arimaCoefficients, maxPQ, Y', E', X', -regressCoefficients)';
%
%  Note:
%
%  To improve runtime performance, the following code segment has been replaced 
%  by a private/undocumented C-MEX file called above, which does exactly the 
%  same thing. To restore the original MATLAB code, simply uncomment the 
%  following code segment:


%    I = ones(numPaths,1);
%
%    for t = (maxPQ + 1):T
%        data   = [I  Y(:,t - LagsAR)  E(:,t - LagsMA)];
%        E(:,t) = data * arimaCoefficients  -  regressCoefficients * X(:,t);
%    end

end

%
% Evaluate the Gaussian negative log-likelihood objective.
%

nLogLikelihoods = 0.5 * (log(2*pi*V(:,(maxPQ + 1):T)) + ((E(:,(maxPQ + 1):T).^2)./V(:,(maxPQ + 1):T)));
nLogL           = sum(nLogLikelihoods, 2);

%
% Trap any degenerate negative log-likelihood values.
%

i = isnan(nLogL) | (imag(nLogL) ~= 0) | isinf(nLogL);

if any(i)
   nLogL(i) = 1e20;
end

end % Gaussian negative log-likelihood function

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function [nLogL, nLogLikelihoods, E] = nLogLikeT(coefficients, Y, X, E, V, LagsAR, LagsMA, numPaths, maxPQ, T)

%
% Infer the residuals.
%

if isempty(X)

   arimaCoefficients = coefficients(1:(end - 1));

   E = internal.econ.arimaxMex(LagsAR, LagsMA, arimaCoefficients, maxPQ, Y', E')';

%
%  Note:
%
%  To improve runtime performance, the following code segment has been replaced 
%  by the private/undocumented C-MEX file called above, which does exactly the 
%  same thing. To restore the original MATLAB code, simply uncomment the 
%  following code segment:


%    I = ones(numPaths,1);
% 
%    for t = (maxPQ + 1):T
%        data   = [I  Y(:,t - LagsAR)  E(:,t - LagsMA)];
%        E(:,t) = data * arimaCoefficients;
%    end

else

   arimaCoefficients   = coefficients(1:(end - size(X,1) - 1));
   regressCoefficients = coefficients((end - size(X,1)):(end - 1))';

   E = internal.econ.arimaxMex(LagsAR, LagsMA, arimaCoefficients, maxPQ, Y', E', X', -regressCoefficients)';

%
%  Note:
%
%  To improve runtime performance, the following code segment has been replaced 
%  by the private/undocumented C-MEX file called above, which does exactly the 
%  same thing. To restore the original MATLAB code, simply uncomment the 
%  following code segment:


%    I = ones(numPaths,1);
% 
%    for t = (maxPQ + 1):T
%        data   = [I  Y(:,t - LagsAR)  E(:,t - LagsMA)];
%        E(:,t) = data * arimaCoefficients  -  regressCoefficients * X(:,t);
%    end

end

%
% Evaluate standardized Student's t negative log-likelihood function.
%

DoF = coefficients(end);

if DoF <= 200

% 
%  Standardized t-distributed residuals.
%

   nLogLikelihoods =  0.5 * (log(V(:,(maxPQ + 1):T))  +  (DoF + 1) * log(1 + (E(:,(maxPQ + 1):T).^2)./(V(:,(maxPQ + 1):T) * (DoF - 2))));
   nLogLikelihoods =  nLogLikelihoods - log(gamma((DoF + 1)/2) / (gamma(DoF/2) * sqrt(pi * (DoF - 2))));
   nLogL           =  sum(nLogLikelihoods, 2);

else

% 
%  Gaussian residuals.
%

   nLogLikelihoods = 0.5 * (log(2*pi*V(:,(maxPQ + 1):T)) + ((E(:,(maxPQ + 1):T).^2)./V(:,(maxPQ + 1):T)));
   nLogL           = sum(nLogLikelihoods, 2);

end

%
% Trap any degenerate negative log-likelihood values.
%

i = isnan(nLogL) | (imag(nLogL) ~= 0) | isinf(nLogL);

if any(i)
   nLogL(i) = 1e20;
end

end % t-distribution negative log-likelihood function

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function Mdl = loadobj(OBJ)
 
%
% The following logic relies on the fact that the "Description" property was
% added to ARIMA models in MATLAB release R2018a.
%

if any(strcmp('Description', fieldnames(OBJ)))
%
%  Since the input structure (OBJ) has a "Description" field, the model was 
%  saved in a recent release and so we simply assign the input structure (OBJ)
%  to the output model (Mdl).
%
   Mdl = OBJ;
else
%
%  Since the input structure (OBJ) does NOT have a "Description" field, the
%  model was saved in a release prior to R2018a and we must create the output 
%  model (Mdl) from the components of the input structure (OBJ).
%
   Mdl              = arima();
   Mdl.Constant     = OBJ.Constant;
   Mdl.LagOpLHS{1}  = OBJ.LagOpLHS{1};
   Mdl.LagOpLHS{2}  = OBJ.LagOpLHS{2};
   Mdl.LagOpRHS{1}  = OBJ.LagOpRHS{1};
   Mdl.LagOpRHS{2}  = OBJ.LagOpRHS{2};
   Mdl.D            = OBJ.D;
   Mdl.Seasonality  = OBJ.Seasonality;
   Mdl.Beta         = OBJ.Beta;
   Mdl.Variance     = OBJ.Variance;
   Mdl.Distribution = OBJ.Distribution;
   Mdl.Description  = getModelSummary(Mdl);
end
 
end % Load Object

end % Methods (Static, Hidden)

end % Class definition