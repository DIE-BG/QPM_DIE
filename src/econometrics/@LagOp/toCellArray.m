function [coefficients,lags] = toCellArray(A)
%TOCELLARRAY Convert a lag operator polynomial object to a cell array
%
% Syntax:
%
%	[coefficients,lags] = toCellArray(A)
%
% Description:
%
%   Convert a lag operator polynomial object A(L) to an equivalent cell
%   array.
%
% Inputs Arguments:
%
%	A - Lag operator polynomial object A(L), as produced by LagOp.
%
% Output Arguments:
%
%   coefficients - Cell array equivalent to the lag operator polynomial
%     A(L). 
%
%   lags - Vector of unique integer lags associated with the polynomial 
%     coefficients. Elements of lags are in ascending order. The first
%     element of lags is the smaller of the smallest nonzero coefficient
%     lag of the object and zero; the last element of lags is the degree of
%     the polynomial. That is, lags = [min(A.Lags,0), 1, 2, ... A.Degree].
%
% Notes:
%
%   o LagOp objects implicitly store polynomial lags and corresponding 
%     coefficient matrices of zero-valued coefficients via lag-based
%     indexing. However, cell arrays conform to traditional element
%     indexing rules, and must explicitly store zero coefficient matrices.
%
%   o The output cell array is equivalent to the input lag operator
%     polynomial in the sense that the same lag operator is created when
%     the output coefficients and lags are used to create a new LagOp
%     object. That is, the following two statements produce the same
%     polynomial A(L):
%
%       [coefficients,lags] = toCellArray(A);
%        A = LagOp(coefficients,'Lags',lags);
%
% See also LagOp/LagOp.

% Copyright 2018 The MathWorks, Inc.

% Perform the polynomial-to-cell-array conversion:

lags         = union(A.Lags, 0);  % Ensure zero is included
lags         = min(lags):max(lags);
dimension    = A.Dimension;
coefficients = mat2cell([A.Coefficients{lags(:)}], dimension, repmat(dimension,1,numel(lags)));

end