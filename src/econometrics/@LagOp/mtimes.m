function C = mtimes(A,B,varargin)
%MTIMES Lag operator polynomial multiplication
%
% Syntax:
%
%   C = A*B
%   C = mtimes(A,B,'Tolerance',tolerance)
%
% Description:
%
%   Given two lag operator polynomials, A(L) and B(L), perform a polynomial
%   multiplication so that C(L) = A(L)*B(L).
%
% Input Arguments:
%
%   A - First lag operator polynomial object factor, as produced by LagOp,
%     in the product A(L)*B(L).
%
%   B - Second lag operator polynomial object factor, as produced by LagOp,
%     in the product A(L)*B(L).
%
%   If at least one of A or B is a lag operator polynomial object, the
%   other can be a cell array of matrices (polynomial coefficients at lags 
%   0, 1, 2, ...), or a single matrix (zero-degree polynomial).
%
% Optional Input Parameter Name/Value Pairs:
%
%   'Tolerance'  Nonnegative scalar tolerance used to determine which
%                coefficients are included in the result. The default
%                tolerance is 1e-12. Specifying a tolerance greater than 0
%                allows the user to exclude polynomial lags with near-zero
%                coefficients. A coefficient matrix of a given lag is
%                excluded only if the magnitudes of all elements of the
%                matrix are less than or equal to the specified tolerance.
%
% Output Argument:
%
%   C - Product lag operator polynomial object, with C(L) = A(L)*B(L). 
%
% Notes:
%
%   o The multiplication operator (*) invokes MTIMES, but the optional
%     coefficient tolerance is available only by calling MTIMES directly. 
%
% See also LagOp/LagOp, LagOp/mrdivide, LagOp/mldivide. 

% Copyright 2010 The MathWorks, Inc.

% Check input arguments:

if nargin > 4
    
   error(message('econ:LagOp:mtimes:TooManyInputs'))
end

parser = inputParser;
parser.addRequired('A');
parser.addRequired('B');
parser.addParamValue('Tolerance',LagOp.ZeroTolerance,@LagOp.checkTolerance);
parser.parse(A,B,varargin{:});

A         = parser.Results.A;
B         = parser.Results.B;
tolerance = parser.Results.Tolerance;

if isempty(A) || isempty(B)
    
   error(message('econ:LagOp:mtimes:MissingInputs'))
end

% Allow either input to be a matrix, a cell array, or some other object, 
% and attempt a LagOp conversion. The tolerance is not applied here.

if ~isa(A,'LagOp')
   A = LagOp(A,'Tolerance',0); % Convert A(L) to a LagOp object
elseif ~isa(B,'LagOp')
   B = LagOp(B,'Tolerance',0); % Convert B(L) to a LagOp object
end

if A.Dimension ~= B.Dimension
    
   error(message('econ:LagOp:mtimes:InconsistentDimension'))
     
end

% Determine the nonzero coefficient lags found in the polynomial product:
%
%  o The call to MESHGRID converts any lags associated with nonzero
%    coefficients of A(L) and B(L) to compatible 2D "grid matrices" of the
%    same size. 
%
%  o These matrices are then converted to column vectors in which each
%    element of the A(L) polynomial vector forms a lag-indexed pair with
%    the corresponding element of the B(L) polynomial vector.
%
%  o The row-wise sum of these A(L) + B(L) polynomial vectors then
%    indicates the lag of the product polynomial vector C(L) to which the
%    corresponding coefficient matrix product contributes. In general,
%    there will be more than one matrix product associated with any given
%    lag of the product polynomial, C(L).

[A_Lags,B_Lags] = meshgrid(A.Lags,B.Lags); % Convert to ND grids

A_Lags = A_Lags(:);              % Lags found in A(L)
B_Lags = B_Lags(:);              % Lags found in B(L)
C_Lags = A_Lags + B_Lags;        % Lags found in C(L)
lags   = unique([0 ; C_Lags])';  % Retain just the unique lags, and include 0

% Now determine the actual coefficients of the product polynomial:
%
% o The outer loop iterates over all product lags (each lag is guaranteed
%   to be in the product).
%
% o The inner loop accumulates all individual coefficient products of
%   A(L)*B(L) which contribute to a given lag of C(L).

coefficients(1:numel(lags)) = {zeros(A.Dimension)};

for L = 1:numel(lags)

    rows = find(lags(L) == C_Lags); % Identify rows which contribute to C(L)

    for iRow = 1:numel(rows)       % Accumulate the sum
        lagA = A_Lags(rows(iRow)); % Lag of A(L) included in this term
        lagB = B_Lags(rows(iRow)); % Lag of B(L) included in this term
        product = A.Coefficients{lagA} * B.Coefficients{lagB};
        coefficients{L} = coefficients{L} +  product;
    end

end

% Create the LagOp object and apply the tolerance:

C = LagOp(coefficients,'Lags',lags,'Tolerance',tolerance);

end