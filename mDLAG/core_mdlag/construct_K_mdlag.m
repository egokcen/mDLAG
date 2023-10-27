function [K_big] = construct_K_mdlag(params,T)
%
% [K_big] = construct_K_mdlag(params, T)
%
% Description: Constructs full GP covariance matrix across all latent
%              state dimensions and timesteps.
%
% Arguments:
%
%     params  -- Structure containing mDLAG model parameters.
%                Contains the following relevant fields:
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma   -- (1 x xDim) array; GP timescales in units of 
%                               time are given by 'binWidth ./ sqrt(gamma)'                                                    
%                    eps     -- (1 x xDim) array; GP noise variances
%                    D       -- (numGroups x xDim) array; delays from 
%                               latents to observed variables. NOTE: Delays
%                               are reported as (real-valued) number of 
%                               time-steps.
%                    xDim    -- int; number of latent variables
%                    yDims   -- (1 x numGroups) array; dimensionalities of 
%                               each observed group
%
%     T       -- int; number of timesteps
%
% Outputs:
%
%     K_big   -- (xDim * numGroups * T) x (xDim * numGroups * T) array;
%                GP covariance matrix            
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     27 Sep 2022 -- Initial full revision.

xDim         = params.xDim;
yDims        = params.yDims;
numGroups    = length(yDims);
K_big        = zeros(xDim*numGroups*T);
mT           = numGroups*T;

Tdif = repmat(1:T,numGroups,1); % (m x T)
Tdif = repmat(Tdif(:)',mT,1) - repmat(Tdif(:),1,mT); % (mT x mT)
for j = 1:xDim
    Delayall = params.D(:,j); % (m x 1)
    Delaydif = repmat(Delayall,T,1);    % (mT x 1)
    Delaydif = repmat(Delaydif',mT,1) - repmat(Delaydif,1,mT); % (mT x mT)
    deltaT = Tdif - Delaydif; 
    deltaTsq = deltaT.^2;
    switch(params.covType)
        case 'rbf'
            temp = exp(-0.5*params.gamma(j)*deltaTsq);          
    end
    K_j = (1-params.eps(j))*temp + params.eps(j)*eye(mT); % (mT x mT)
    
    idx = j:xDim:xDim*numGroups*T;
    K_big(idx,idx) = K_j;
end
