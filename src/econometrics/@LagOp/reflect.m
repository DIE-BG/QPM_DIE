function B = reflect(A)
%REFLECT Reflect lag operator polynomial coefficients around lag zero
%
% Syntax:
%
%   B = reflect(A)
%
% Description:
%
%   Given a lag operator polynomial object A(L), negate all coefficient
%   matrices except the coefficient matrix at lag 0. For example, given a
%   polynomial of degree p,
%
%     A(L) = A0 + A1*L + A2*L^2 + ... + AP*L^p
%
%   the reflected polynomial B(L) is
%
%     B(L) = A0 - A1*L - A2*L^2 - ... - AP*L^p
%
%   with the same degree and dimension as A(L).
%
% Input Arguments:
%
%   A - Lag operator polynomial object A(L), as produced by LagOp.
%
% Output Arguments:
%
%   B - Lag operator polynomial B(L), the reflection of A(L). 
%
% See also LagOp/LagOp.

% Copyright 2018 The MathWorks, Inc.

% Reflect the polynomial A(L):

lags = union(A.Lags, 0);
A = toCellArray(A);             % Convert to cell vector for indexing performance (avoids SUBSREF)

for L = 1:numel(lags)
    if lags(L) > 0
       A{lags(L) + 1} = -A{lags(L) + 1};
    end
end

% Create the reflected LagOp object:

B = LagOp(A, 'Tolerance',0);

end