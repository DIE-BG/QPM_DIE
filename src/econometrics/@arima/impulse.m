function varargout = impulse(Mdl, varargin)
%IMPULSE Impulse response (dynamic multipliers) of an ARIMA model
%
% Syntax:
%
%   Y = impulse(Mdl)
%   Y = impulse(Mdl,numObs)
%   impulse(...)
%
% Description:
%
%   Given a univariate ARIMA model, compute the dynamic response of the
%   system to a single impulse, or innovation shock, of unit size. The 
%   specific impulse response calculated here is referred to as the dynamic 
%   multiplier, defined as the partial derivative of the output response
%   with respect to an innovation shock at time zero. Expressed as a function
%   of time, the sequence of dynamic multipliers measures the sensitivity 
%   of the output to a purely transitory change in the innovations process.
%   If no output is requested, a discrete stem plot of the impulse response 
%   is created in the current figure window.
%
% Input Arguments:
%
%   Mdl - ARIMA model specification object, as produced by the ARIMA 
%     constructor or ARIMA/ESTIMATE method.
%
% Optional Inputs:
%
%   numObs - Positive integer indicating the number of observations included
%     in the impulse response (the number of periods for which the impulse 
%     response is computed). If unspecified, the number of observations is
%     determined by the polynomial division algorithm of the underlying Lag 
%     Operator Polynomials (see LagOp/mldivide for details).
%
% Output Arguments:
%
%   Y - numObs-by-1 column vector of the impulse response of the model Mdl.
%
% Notes:
%
%   o For observed model responses y(t) and innovations e(t), the output 
%     response y(t) is computed by shocking the system with a unit impulse 
%     such that e(0) = 1, with all past observations of y(t) and all past 
%     and future shocks of e(t) set to zero. This allows the examination of 
%     the dynamic behavior of the system to a one-time, unit impulse of the 
%     input innovations e(t). Additionally, since the response is viewed as 
%     the partial derivative with respect to an innovation shock at time 
%     zero, the presence of a model constant has no effect on the output.
%
%     Since the innovations e(t) may be interpreted as the 1-step-ahead 
%     forecast errors, the specific impulse response computed here is also 
%     sometimes referred to as the "forecast error impulse response".
%
%   o When the number of observations (see numObs above) is specified, the 
%     impulse response is computed by filtering a unit impulse followed by
%     a vector of zeros of appropriate length; the filtering algorithm is 
%     very fast and results in an impulse response of known length. 
%
%     However, when numObs is unspecified the model is converted to a 
%     truncated infinite-degree moving average representation using the 
%     relatively slow Lag Operator Polynomial division algorithm and results 
%     in an impulse response of generally unknown length.
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
%   [4] Lutkepohl, H. New Introduction to Multiple Time Series Analysis. 
%       New York, NY: Springer-Verlag, 2007.
%
% See also ARIMA, SIMULATE, FILTER.

% Copyright 2018 The MathWorks, Inc.

if numel(Mdl) > 1
   error(message('econ:arima:impulse:NonScalarModel'))
end

%
% Check input parameters and set defaults.
%

parser = inputParser;
parser.addOptional('numObs', NaN, @(x) validateattributes(x, {'double'}, {'scalar' 'integer' '>' 0}, '', 'number of observations'));

try
  parser.parse(varargin{:});
catch exception
  exception.throwAsCaller();
end

numObs = parser.Results.numObs;

%
% Determine which method is used to compute the impulse response. 
%

wState  = warning;                           % Save warning state
cleanUp = onCleanup(@() warning(wState));    % Restore warning state

if isnan(numObs)  % The number of observations specified?

%
%  The user did not specify the number of observations, so convert the model 
%  to a truncated infinite-degree MA representation by performing a left 
%  division of Lag Operator Polynomials. This approach is relatively slow 
%  and results in an impulse response of (generally) unknown length.
%

   AR = Mdl.getLagOp('Compound AR');

   if any(isnan(cell2mat(toCellArray(AR))))
      error(message('econ:arima:impulse:UnspecifiedAR'))
   end

   MA = Mdl.getLagOp('Compound MA');

   if any(isnan(cell2mat(toCellArray(MA))))
      error(message('econ:arima:impulse:UnspecifiedMA'))
   end
   
   warning('off', 'econ:LagOp:mldivide:WindowNotOpen')   
   warning('off', 'econ:LagOp:mldivide:WindowIncomplete')

   try
     Y = cell2mat(toCellArray(AR \ MA));        % Infinite-degree MA representation
   catch exception
     exception.throwAsCaller();
   end

else

%
%  The user specified the number of observations, so filter a unit impulse 
%  assuming all pre-sample data are zero. Since the FILTER method assumes 
%  the input noise process is a standardized process, the output is then 
%  normalized by the first response to ensure that Y(0) = 1. This direct
%  filtering algorithm is very fast and results in an impulse response of 
%  known length. 
%

   try
     Mdl.PrivateVariance = 1;    % So no subsequent scaling is required
     Mdl.PrivateConstant = 0;    % So any constant is effectively ignored
     Y            = filter(Mdl, [1 ; zeros(numObs - 1,1)], 'Y0', zeros(Mdl.P,1), 'Z0', zeros(Mdl.Q,1));
   catch exception
     exception.throwAsCaller();
   end

end

%
% Plot the impulse response if necessary.
%

if nargout == 0
   stem((0:numel(Y)-1)', Y, 'filled');
   xlabel('Observation Time')
   ylabel('Response')
   title('Impulse Response')
else
   varargout = {Y(:)};
end

end