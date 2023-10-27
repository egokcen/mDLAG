function Kx = construct_Kx_mdlag(params,T)
%
% Kx = construct_Kx_mdlag(params, T)
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
%     Kx -- (xDim x 1) cell array; Kx{j} is a (numGroups*T) x (numGroups*T)
%           GP covariance matrix for latent j           
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     01 Feb 2023 -- Initial full revision.

xDim         = params.xDim;
yDims        = params.yDims;
numGroups    = length(yDims);
mT           = numGroups*T;
Kx = cell(xDim,1);
for j = 1:xDim
    Kx{j} = zeros(mT); % M -> T
end

Tdif = repmat(1:T,numGroups,1)'; % (T x m)
Tdif = repmat(Tdif(:)',mT,1) - repmat(Tdif(:),1,mT); % (mT x mT), M -> T
for j = 1:xDim
    Delaydif = repmat(params.D(:,j),1,T)'; % (T x m)
    Delaydif = repmat(Delaydif(:)',mT,1) - repmat(Delaydif(:),1,mT); % (mT x mT), M -> T
    deltaT = Tdif - Delaydif; % M -> T
    deltaTsq = deltaT.^2;
    switch(params.covType)
        case 'rbf'
            temp = exp(-0.5*params.gamma(j)*deltaTsq);          
    end
    Kx{j} = (1-params.eps(j))*temp + params.eps(j)*eye(mT); % (mT x mT), % M -> T
end
