function C = minus(A,B,varargin)
%MINUS Lag operator polynomial subtraction
%
% Syntax:
%
%   C = minus(A,B,'Tolerance',tolerance)
%   C = A - B
%
% Description:
%
%   Given two lag operator polynomials, A(L) and B(L), perform a polynomial
%   subtraction so that C(L) = A(L) - B(L).
%
% Inputs Arguments:
%
%   A - First (minuend) lag operator polynomial object, as produced by
%     LagOp, in the difference A(L) - B(L).
%
%   B - Second (subtrahend) lag operator polynomial object, as produced by
%     LagOp, in the difference A(L) - B(L).
%
%   If at least one of A or B is a lag operator polynomial object, the
%   other can be a cell array of matrices (polynomial coefficients at lags
%   0, 1, 2, ...), or a single matrix (zero-degree polynomial). 
%
% Optional Input Parameter Name/Value Pairs:
%
%   'Tolerance'  Nonnegative scalar tolerance used to determine which lags 
%                are included in the difference. The default is 1e-12.
%                Specifying a tolerance greater than the default excludes
%                lags with near-zero coefficients. A lag is excluded if the
%                magnitudes of all elements of the coefficient matrix are
%                less than or equal to the specified tolerance.
%
% Output Argument:
%
%   C - Difference lag operator polynomial object, with C(L) = A(L)-B(L). 
%
% Note:
%
%   o The subtraction operator (-) invokes MINUS, but the optional
%     coefficient tolerance is available only by calling MINUS directly. 
%
% See also LagOp/LagOp, LagOp/plus.

% Copyright 2010 The MathWorks, Inc. 

% Check input arguments:

if nargin > 4
    
   error(message('econ:LagOp:minus:TooManyInputs'))
     
end

parser = inputParser;
parser.addRequired('A');
parser.addRequired('B');
parser.addParamValue('Tolerance', LagOp.ZeroTolerance,@LagOp.checkTolerance);
parser.parse(A,B,varargin{:});

A         = parser.Results.A;
B         = parser.Results.B;
tolerance = parser.Results.Tolerance;

if isempty(A) || isempty(B)
    
   error(message('econ:LagOp:minus:MissingInputs'))
     
end

% Allow either input to be a matrix, a cell array, or some other object, 
% and attempt a LagOp conversion. The tolerance is not applied here.

if ~isa(A,'LagOp')
   A = LagOp(A,'Tolerance',0); % Convert A(L) to a LagOp object
elseif ~isa(B,'LagOp')
   B = LagOp(B,'Tolerance',0); % Convert B(L) to a LagOp object
end

if A.Dimension ~= B.Dimension
    
   error(message('econ:LagOp:minus:InconsistentDimension'))
     
end

% Perform the polynomial subtraction:

lags = union([0 A.Lags],B.Lags);
coefficients = cell(numel(lags),1);

for L = 1:numel(lags)
    coefficients{L} = A.Coefficients{lags(L)} - B.Coefficients{lags(L)};
end

% Create the LagOp object and apply the tolerance:

C = LagOp(coefficients,'Lags',lags,'Tolerance',tolerance);

end