function [m0, l0, a0, b0] = DCM_get_prior_all(DCM, options)
% [m0, l0, a0, b0] = DCM_get_prior_all(DCM)
% 
% Returns prior parameters on model parameters (theta) and noise precision
% (tau) for the full connectivity model. Necessary for the sparse version
% of rDCM.
% 
%   Input:
%   	DCM             - model structure
%
%   Output:
%       m0              - prior mean
%       l0              - prior covariance
%       a0              - prior shape parameter
%       b0              - prior rate parameter
%
 
% ----------------------------------------------------------------------
% 
%   Authors: Haifeng Wu (whf5469@gmail.com)
% 
%   Copyright (C) 2023 School of Electricical & Information Technology
%                         Yunnan Minzu University
%
% ----------------------------------------------------------------------


% get the prior mean and covariance from SPM
[pE,pC] = DCM_fmri_priors(ones(size(DCM.a)),ones(size(DCM.b)),ones(size(DCM.cb)),ones(size(DCM.d)));

% number of regions and inputs
[nr, nu] = size(pE.C);

% set the prior mean of endogenous parameters to zero
pE.A = zeros(size(pE.A))+diag(diag(pE.A));

% prior precision
pC.A       = 1./pC.A;
pC.B       = 1./pC.B;
pC.C       = 1./pC.C;
pC.D       = 1./pC.D;
pC.transit = 1./pC.transit;
pC.decay   = 1./pC.decay;
pC.epsilon = 1./pC.epsilon;

% prior precision
l0 = [pC.A pC.C];

% Setting priors on noise
a0 = 2;
b0 = 1;

% prior mean 
if options.type == 'S'
    m0 = [pE.A pE.C];
elseif options.type == 'E'
    m0 = options.m0_coef*[diag(ones(length(DCM.a),1)) DCM.cb];
    %m0 = options.m0_coef*[diag(ones(length(DCM.a),1)) pE.C];
end

end
