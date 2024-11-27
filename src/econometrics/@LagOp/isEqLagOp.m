function indicator = isEqLagOp(A,B,varargin)
%ISEQLAGOP Determine if two LagOp objects are the same mathematical polynomial
%
% Syntax:
%
%   indicator = isEqLagOp(A,B)
%   indicator = isEqLagOp(A,B,'Tolerance',tolerance)
%
% Description:
%
%   Two lag operator polynomials, A(L) and B(L), are equal if they have
%   identical nonzero coefficient matrices at all lags. Equality is
%   evaluated within a specified or default positive tolerance.
%
% Input Arguments:
%
%   A - Lag operator polynomial object, as produced by LagOp, against which
%	    the equality of B is tested.
%
%   B - Lag operator polynomial object, as produced by LagOp, against which
%	    the equality of A is tested.
%
%   If at least one of A or B is a lag operator polynomial object, the
%   other can be a cell array of matrices (polynomial coefficients at lags 
%   0, 1, 2, ...), or a single matrix (zero-degree polynomial).
%
% Optional Input Parameter Name/Value Pair:
%
%   'Tolerance' - Nonnegative scalar tolerance used for testing equality. The
%                 default is 1e-12. Specifying a tolerance greater than the
%                 default relaxes the comparison criterion. Two polynomials
%                 are deemed sufficiently close to indicate equality if the
%                 differences in magnitude of all elements of all coefficient 
%                 matrices at all lags are less than or equal to the specified
%                 tolerance.
%
% Output Argument:
%
%   indicator - Boolean indicator for the equality test. TRUE indicates the
%     two polynomials are identical to within tolerance; FALSE indicates
%     the two polynomials are not identical to within tolerance.
%
% See also LagOp/LagOp.

% Copyright 2010 The MathWorks, Inc.

% Check input arguments:

if nargin > 4
    
    error(message('econ:LagOp:isEqLagOp:TooManyInputs'))
      
end

parser = inputParser;
parser.addRequired('A');
parser.addRequired('B');
parser.addParamValue('Tolerance',LagOp.ZeroTolerance,@LagOp.checkTolerance);
parser.parse(A,B,varargin{:});

A = parser.Results.A;
B = parser.Results.B;
tolerance = parser.Results.Tolerance;

if isempty(A) || isempty(B)
    
   error(message('econ:LagOp:isEqLagOp:MissingInputs'))
     
end

% Allow either input to be a matrix, a cell array, or some other object, 
% and attempt a LagOp conversion. The tolerance is not applied here.

if ~isa(A,'LagOp')
   A = LagOp(A,'Tolerance',0); % Convert A(L) to a LagOp object
elseif ~isa(B,'LagOp')
   B = LagOp(B,'Tolerance',0); % Convert B(L) to a LagOp object
end

if A.Dimension ~= B.Dimension % Are they the same dimension?
   indicator = false;
   return
end

% Test the coefficients for equality to within tolerance:

Lags = union(A.Lags,B.Lags); % Union of lags allows for tolerance issues

C1 = cell(numel(Lags),1);
C2 = C1;

for L = 1:numel(Lags)
    C1{L} = A.Coefficients{Lags(L)};
    C2{L} = B.Coefficients{Lags(L)};
end

indicator = all(all(abs([C1{:}]-[C2{:}]) <= tolerance));

end