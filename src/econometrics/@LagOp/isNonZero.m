function indicator = isNonZero(OBJ,testLags)
%ISNONZERO Find lags associated with nonzero coefficients of LagOp objects
%
% Syntax:
%
%   indicator = isNonZero(A,testLags)
%
% Description:
%
%   Given a vector of candidate lags to test, determine which lags are
%   associated with nonzero coefficients of a lag operator polynomial A(L). 
%
% Inputs Arguments:
%
%	A - Lag operator polynomial object, as produced by LagOp.
%
%   testLags - Vector of unique integer lags to test in A.
%
% Output Argument:
%
%   indicator - Boolean vector of lag inclusion tests, the same length as
%     testLags. TRUE indicates the corresponding test lag is included in A;
%     FALSE indicates the corresponding test lag is not included in A.
%
% See also LagOp/LagOp.

% Copyright 2010 The MathWorks, Inc.    $Date:
% 2009/09/16 08:53:28 $

indicator = false(size(testLags)); % TRUE indicates a lag is specified
currentLags = OBJ.Lags;

for L = 1:numel(testLags)
    indicator(L) = any(currentLags == testLags(L));
end

end