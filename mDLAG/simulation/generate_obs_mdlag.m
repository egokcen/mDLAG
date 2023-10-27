function seq = generate_obs_mdlag(seq, params, varargin)
% 
% seq = generate_obs_mdlag(seq, params,...)
%
% Description: Generate observations given a set of latent sequences, all
%              from the mDLAG generative model.
%
% Arguments:
%
%     Required:
%
%     seq -- structure whose nth entry (corresponding to the nth 
%            sequence) has fields
%                trialId   -- int; unique trial (sequence) identifier  
%                T -- int; number of timesteps
%                (latentfield) -- (numGroups*xDim x T) array; delayed
%                                 latent sequences
%
%     params  -- Structure containing mDLAG model parameters.
%                Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    Cs      -- (1 x numGroups) cell array; List of factor 
%                               loadings 
%                               {(y1Dim x xDim), (y2Dim x xDim), ...}
%                    alphas  -- (numGroups x xDim) array; alpha parameter 
%                               values
%                    phis    -- (1 x numGroups) cell array; List of 
%                               observation precision parameters 
%                               {(y1Dim x 1), (y2Dim x 1), ...}
%                    ds      -- (1 x numGroups) cell array; List of data 
%                               means
%                               {(y1Dim x 1), (y2Dim x 1), ...}
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
%     Optional:
%
%     latentfield -- string; Name of latent data field in seq 
%                    (default: 'xsm')
%     obsfield    -- string; Name of observation data field in seq 
%                    (default: 'y')
%
% Outputs:
%     seq -- structure whose nth entry (corresponding to the nth 
%            sequence) has fields
%                trialId   -- int; unique trial (sequence) identifier  
%                T -- int; number of timesteps
%                (latentfield) -- (numGroups*xDim x T) array; 
%                                 delayed latent sequences
%                (obsfield) -- (yDim x T) array; observation sequence
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     27 Sep 2022 -- Initial full revision.

latentfield = 'xsm';
obsfield = 'y';
extraOpts = assignopts(who, varargin);

N = length(seq);
numGroups = length(params.yDims);
yDim = sum(params.yDims);
obs_block_idxs = get_block_idxs(params.yDims);
lat_block_idxs = get_block_idxs(params.xDim.*ones(1,numGroups));

for n = 1:N
    T = seq(n).T;
    X = seq(n).(latentfield);
    seq(n).(obsfield) = nan(yDim,T);
    for groupIdx = 1:numGroups
        obsBlock = obs_block_idxs{groupIdx};
        latBlock = lat_block_idxs{groupIdx};
        ns = mvnrnd(zeros(1,params.yDims(groupIdx)), diag(params.phis{groupIdx}.^(-1)), T)';
        seq(n).(obsfield)(obsBlock(1):obsBlock(2),:) ...
            = params.Cs{groupIdx} * X(latBlock(1):latBlock(2),:) + repmat(params.ds{groupIdx},1,T) + ns;
    end
end
