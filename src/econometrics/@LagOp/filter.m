function [Y,times] = filter(A,X,varargin)
%FILTER Apply a lag operator polynomial to filter a time series
%
% Syntax:
%
%	[Y,times] = filter(A,X)
%   [Y,times] = filter(A,X,'Initial',X0)
%
% Description:
%
%   Given a lag operator polynomial A(L), apply A(L) to time series data
%   X(t). This is equivalent to applying a linear filter to X(t), producing
%   the filtered output series Y(t) = A(L)X(t).
%
% Input Arguments:
%
%	A - Lag operator polynomial object, as produced by LagOp.
%
%   X - numObs-by-numDims matrix of time series data to which the lag
%     operator polynomial A is applied. The last observation is assumed to
%     be the most recent. numDims is the dimension of A, unless X is a row
%     vector, in which case X is treated as a univariate series. For
%     univariate X, the orientation of the output Y is determined by the
%     orientation of the input X.
%
% Optional Input Parameter Name/Value Pair:
%
%	'Initial'  Presample values of the input time series X(t). If 'Initial'
%	           is unspecified, or if the number of presample values is
%	           insufficient to initialize filtering, values are taken from
%	           the beginning of X, reducing the effective sample size of
%	           the output Y. For convenience, scalar presample values are
%	           expanded to provide all numPresampleObs-by-numDims presample
%	           values, and data is not taken from X. If more presample
%	           values are specified than necessary, only the most recent
%	           values are used. For univariate X, presample values can be a
%	           row or a column vector.
%
% Output Arguments:
%
%	Y - Filtered input series X(t), with Y(t) = A(L)X(t). 
%
%   times - A vector of relative time indices, the same length as Y. Times
%     are expressed relative to, or as an offset from, observations times
%     0, 1, 2, ..., numObs-1 for the input series X(t). For a polynomial of
%     degree p, Y(0) is a linear combination of X(t) for times t = 0, -1,
%     -2, ..., -p (presample data). Y(t) for t > 0 is a linear combination
%     of X(t) for times t, t-1, t-2, ..., t-p.
%
% Note:
%
%   o Filtering is limited to single paths, so matrix data are assumed to
%     be a single path of a multidimensional process, and 3D data (multiple
%     paths of a multidimensional process) are not allowed.
%
% See also LagOp/LagOp.

% Copyright 2010 The MathWorks, Inc.    $Date:
% 2009/09/21 15:15:05 $

if nargin > 4
    
   error(message('econ:LagOp:filter:TooManyInputs'))
     
end

isRow = false; % Default univariate output orientation

if (size(X,1) == 1) && isvector(X)
   X = X(:);     % Guarantee column orientation
   isRow = true; % Set output orientation
end

% Extract some information for input argument consistency checking:

p = A.Degree;    % Degree of lag operator polynomial
n = A.Dimension; % Dimension of lag operator

% Check input arguments:

parser = inputParser;
parser.addRequired ('A');
parser.addRequired ('X', @checkTimeSeries);
parser.addParamValue('Initial',[],@checkInitialValues);
parser.parse(A,X,varargin{:});

A  = parser.Results.A;
X  = parser.Results.X;
X0 = parser.Results.Initial;

if isempty(A)
    
   error(message('econ:LagOp:filter:EmptyPolynomial'))
     
end

if isempty(X)

   % To be consistent with the FILTER function of core MATLAB, an empty
   % input series produces an empty output series:

   Y = [];  times = [];  return
end

% Initialize presample data, if necessary:

if isempty(X0) % This test relies on the default "[]" passed to inputParser.

   % If the X0 is still empty, then no presample data is specified:
   
   if size(X,1) <= p
       
      error(message('econ:LagOp:filter:InsufficientData'))
        
   end

   % In the absence of explicit presample observations, any required
   % initial values are stripped from the beginning of the input time
   % series X(t), reducing the effective sample size of the output filtered
   % time series Y(t):

   nX0 = 0;            % No presample observations specified
   X0 = X(1:p,:);      % Strip initial observations
   X = X((p+1):end,:); % Reduce the size of input series

elseif numel(X0) == 1 % Scalar presample data specified

   nX0 = p;            % All presample observations specified
   X0 = X0(ones(p,n)); % Scalar expansion

else % Non-scalar presample data specified 

   % In the event that both X(t) and X0(t) are specified as vectors, assume
   % a univariate time series:

   if isvector(X) && isvector(X0)
      X0 = X0(:); % Guarantee column orientation
   end

   % If sufficient presample observations exist, retain only what is
   % necessary; otherwise, use all presample observations provided and
   % strip from the beginning of X(t), as necessary:

   nX0 = size(X0,1);

   if (size(X,1) + nX0) <= p
       
      error(message('econ:LagOp:filter:NotEnoughData'))
        
   end

   X0 = [X0((nX0-min(p,nX0)+1):nX0,:); X(1:(p-nX0),:)];
   X = X(max(p-nX0+1,1):end,:);

end


% Extract the nonzero coefficients from A(L) for runtime performance:

Lags = A.Lags; % Get all lags associated with nonzero coefficients
nLags = numel(Lags);
C = cell(1,nLags);

for i = 1:nLags
    C{i} = A.Coefficients{Lags(i)};
end

% Filter the time series: Y(t) = A(L)X(t).
%
% Although the input time series X(t) is a numObs-by-numDims matrix, for
% runtime performance the filtered output series Y(t) is temporarily
% formatted so that time runs across columns rather than down rows. This
% allows a vector observation to be sequentially accessed and stored one
% column at a time.

T = size(X,1);
Y = zeros(n,T);

for t = 1:T
    for i = 1:nLags
        if t <= Lags(i)
           Y(:,t) = Y(:,t) + C{i}*X0(t-Lags(i)+ p,:)';
        else
           Y(:,t) = Y(:,t) + C{i}*X(t-Lags(i),:)';
        end
    end
end

if ~isRow
   Y = Y.'; % Transpose to column-oriented format
end

if nargout > 1 % Are relative times requested?
   times = (max(p-nX0,0):max(p-nX0,0)+(T-1))';
   if isRow
      times = times';
   end
end

%-------------------------------------------------------------------------
% Check input X
function OK = checkTimeSeries(X)

   if isnumeric(X)
       
      if (size(X,2) ~= n)
          
         error(message('econ:LagOp:filter:checkTimeSeries:IncompatibleDimension'))
           
      end
      
      if ndims(X) > 2
          
         error(message('econ:LagOp:filter:checkTimeSeries:Non2DMatrix'))
           
      end
      
   else
       
      error(message('econ:LagOp:filter:checkTimeSeries:InvalidDataTypeX'))
        
   end

   OK = true;

end % checkTimeSeries

%-------------------------------------------------------------------------
% Check input initial values
function OK = checkInitialValues(X0)

   if isnumeric(X0) 

      % In the event that both X(t) and X0(t) are specified as vectors,
      % assume a univariate time series:

      if isvector(X) && isvector(X0)
         X0 = X0(:); % Guarantee column orientation
      end
      
      if (size(X0,2) ~= n) && (numel(X0) ~= 1) % Allow for scalar expansion
          
         error(message('econ:LagOp:filter:checkInitialValues:IncompatibleDimension'))
           
      end
      
      if ndims(X0) > 2
          
         error(message('econ:LagOp:filter:checkInitialValues:Non2DMatrix'))
           
      end
      
   else
       
      error(message('econ:LagOp:filter:checkInitialValues:InvalidDataTypeX0'))
        
   end

   OK = true;

end  % checkInitialValues

end