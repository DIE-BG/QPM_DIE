function varargout = aicbic(varargin)
%AICBIC Information criteria
%
% Syntax:
%
%   aic = aicbic(logL,numParam)
%   [aic,bic] = aicbic(logL,numParam,numObs)
%   [aic,bic] = aicbic(logL,numParam,numObs,'Normalize',true)
%   [aic,bic,ic] = aicbic(logL,numParam,numObs)
%   [aic,bic,ic] = aicbic(logL,numParam,numObs,'Normalize',true)
%
%
% Description:
%
%   Given loglikelihood values logL obtained by fitting a model to data,
%   compute information criteria to assess model adequacy. Information
%   criteria rank models using measures that balance goodness of fit with
%   parameter parsimony. Models with lower criteria values are preferred.
%
% Input Arguments:
%
%   logL - Loglikelihoods associated with parameter estimates of different
%          models, specified as a vector of numeric values.
%
%   numParam - Number of estimated parameters in the models, specified as a
%          positive integer applied to all elements in logL, or a vector of
%          positive integers having the same length as logL.
%
% Optional Input Argument:
%
%   numObs - Sample sizes used in estimation, specified as a positive
%          integer applied to all elements in logL, or a vector of positive
%          integers having the same length as logL. AICBIC requires numObs
%          for all criteria except the Akaike information criterion, or if
%          'Normalize' is true.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME            VALUE
%
%   'Normalize'     Flag to normalize results by numObs, specified as a
%                   logical value. When true, all output arguments are
%                   divided by numObs. The default is false.
%
% Output Arguments:
%
%   aic - Vector of Akaike information criteria corresponding to elements
%        of logL.
%
%   bic - Vector of Bayesian (Schwarz) information criteria corresponding
%        to elements of logL.
%
%   ic - Structure array with fields:
%
%        aic	Akaike information criteria (AIC)
%        bic	Bayesian (Schwarz) information criteria (BIC)
%        aicc   Corrected Akaike information criteria (AICc)
%        caic	Consistent Akaike information criteria (CAIC)
%        hqc	Hannan-Quinn criteria (HQC)
%
%        ic.aic and ic.bic are the same values returned in aic and bic.
%        AICBIC computes unnormalized criteria as follows:
% 
%        o	AIC  = -2*logL + 2*numParam 
% 
%        o	BIC = -2*logL + log(numObs)*numParam 
% 
%        o	AICC = AIC + [2*numParam*(numParam+1)]/(numObs-numParam-1)
% 
%        o	CAIC = -2*logL + (log(numObs)+1)*numParam
% 
%        o	HQC = -2*logL + 2*log(log(numObs))*numParam
% 
% Notes:
%
%   o Misspecification tests LMTEST, LRATIOTEST, and WALDTEST compare the
%     loglikelihoods of two competing nested models. By contrast, AICBIC
%     accepts the loglikelihoods of individual model fits and returns
%     approximate measures of "information loss" with respect to the data-
%     generating process. Information criteria provide relative rankings of
%     any number of competing models, including non-nested models.
%
%   o In small samples, AIC tends to overfit. To address overfitting, AICc
%     adds a size-dependent correction term that increases the penalty on
%     the number of parameters. AICc approaches AIC asymptotically.
%     Analysis in [3] suggests using AICc when numObs/numParam < 40.
%     
%   o When econometricians compare models with different numbers of
%     autoregressive lags or different orders of differencing, they often
%     scale information criteria by the number of observations [5]. To do
%     this, set numObs to the effective sample size of each estimate, and
%     set 'Normalize' to true.
%
% Example:
%
%   % Simulate DGP
%
%   T = 100;
%   DGP = arima('Constant',1,'AR',[0.2,-0.4],'Variance',1);
%   y = simulate(DGP,T);
% 
%   % Competing models
% 
%   Mdl1 = arima('ARLags',1);
%   Mdl2 = arima('ARLags',1:2);
%   Mdl3 = arima('ARLags',1:3);
% 
%   % Compute log-likelihoods
% 
%   logL = zeros(3,1);
%   [~,~,logL(1)] = estimate(Mdl1,y,'Display','off');
%   [~,~,logL(2)] = estimate(Mdl2,y,'Display','off');
%   [~,~,logL(3)] = estimate(Mdl3,y,'Display','off');
% 
%   % Compute and compare information criteria
% 
%   numParam = [3;4;5];
%   numObs = T*ones(3,1);
%   [~,~,ic] = aicbic(logL,numParam,numObs)
% 
% References:
%
%   [1] Akaike, H. "Information Theory and an Extension of the Maximum
%       Likelihood Principle." In: Petrov B., Csaki F., editors. Second
%       International Symposium on Information Theory. Budapest: Akademiai
%       Kiado, 1973, pp. 267-281.
%
%   [2] Akaike, H. "A New Look at the Statistical Model Identification."
%       IEEE Transactions on Automatic Control. Vol. 19, No, 6, 1974, 
%       pp. 716-723.
% 
%   [3] Burnham, K. and D. Anderson. Model Selection and Multimodel
%       Inference: A Practical Information-Theoretic Approach, 2nd Ed.
%       New York: Springer, 2003.
% 
%   [4] Hannan, E. and B. Quinn. "The Determination of the Order of an
%       Autoregression." Journal of the Royal Statistical Society, Series
%       B. Vol. 41, 1979, pp. 190-195.
% 
%   [5] Lutkepohl, H. and M. Kratzig. Applied Time Series Econometrics.
%       Cambridge: Cambridge University Press, 2004.
% 
%   [6] Schwarz, G. "Estimating the Dimension of a Model." The Annals of
%       Statistics. Vol. 6, No. 2, 1978, pp. 461-464.
%
% See also LMTEST, LRATIOTEST, WALDTEST.

% Copyright 2020 The MathWorks, Inc.

% Parse inputs and set defaults:

parseObj = inputParser;

addRequired(parseObj,'logL',...
	@(x)validateattributes(x,{'numeric'},{'vector'}))
    
addRequired(parseObj,'numParam',...
	@(x)validateattributes(x,{'numeric'},{'vector','positive','integer'}))
    
addOptional(parseObj,'numObs',1,...
	@(x)validateattributes(x,{'numeric'},{'vector','positive','integer'}))
   
addParameter(parseObj,'Normalize',false,...
	@(x)validateattributes(x,{'numeric','logical'},{'scalar','binary'}))
   
parse(parseObj,varargin{:});

logL = parseObj.Results.logL;
numParam = parseObj.Results.numParam;
numObs = parseObj.Results.numObs;
normalizeOutput = parseObj.Results.Normalize;

rowOutput = (size(logL,1) == 1); % Output orientation

logL = logL(:)';
numParam = numParam(:)';
numObs = numObs(:)';

if (length(numParam) > 1) && (length(numParam) ~= length(logL))
       
	error(message('econ:aicbic:VectorLengthMismatch1'))
        
elseif length(numParam) == 1
     
	numParam = numParam(ones(1,length(logL))); % Scalar expansion
      
end

if (length(numObs) > 1) && (length(numObs) ~= length(logL))
       
	error(message('econ:aicbic:VectorLengthMismatch2'))
        
elseif length(numObs) == 1
     
	numObs = numObs(ones(1,length(logL))); % Scalar expansion
      
end

aic  = -2*logL + 2*numParam;
bic  = -2*logL + log(numObs).*numParam; 
aicc = aic + (2*numParam.*(numParam +1))./(numObs-numParam-1);
caic = -2*logL + (log(numObs)+1).*numParam;
hqc  = -2*logL + 2*log(log(numObs)).*numParam;

if normalizeOutput
    
    aic  = aic./numObs;
    bic  = bic./numObs;
    aicc = aicc./numObs;
    caic = caic./numObs;
    hqc  = hqc./numObs;    
    
end

ic.aic  = aic;
ic.bic  = bic;
ic.aicc = aicc;
ic.caic = caic;
ic.hqc  = hqc;

if ~rowOutput
   
    aic = aic';
    bic = bic';
    
end

varargout{1} = aic;
varargout{2} = bic;
varargout{3} = ic;