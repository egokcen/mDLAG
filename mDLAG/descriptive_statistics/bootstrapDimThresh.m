function [R2, MSE] = bootstrapDimThresh(seq, params, cutoff_sharedvar, cutoff_snr, numBoot, varargin)
%
% [R2, MSE] = bootstrapDimThresh(seq, params, cutoff_sharedvar, cutoff_snr, numBoot, ...)
%
% Description: Evaluate leave-group-out predictive performance across a
%              specified range of shared variance cutoffs using a
%              non-parametric bootstrap procedure.
% 
% Arguments:
%
%     Required:
%
%     seq     -- data structure, whose nth entry (corresponding to the nth 
%                trial) has fields
%                    trialId      -- unique trial identifier
%                    T (1 x 1)    -- number of timesteps
%                    y (yDim x T) -- observed data
%     params  -- Structure containing mDLAG model parameters.
%                Contains the fields
%         covType    -- string; type of GP covariance (e.g., 'rbf')
%         gamma      -- (1 x xDim) array; GP timescales in units of 
%                       time are given by 'binWidth ./ sqrt(gamma)'                                                    
%         eps        -- (1 x xDim) array; GP noise variances
%         D          -- (numGroups x xDim) array; delays from latents to 
%                       observed variables. NOTE: Delays are reported as 
%                       (real-valued) number of time-steps.
%         d.mean     -- (yDim x 1) array; posterior mean of mean parameter
%         d.cov      -- (yDim x 1) array; diagonal elements of the
%                       posterior covariance matrix of d
%         C.means    -- (numGroups x 1) cell array; yDims(groupIdx) x xDim
%                       mean loadings matrix for each group
%         C.covs     -- (numGroups x 1) cell array; C.covs{groupIdx) is a
%                       (yDims(groupIdx) x 1) cell array, and each element
%                       is a (xDim x xDim) matrix giving the posterior
%                       covariance of a row of C.
%         C.moments  -- (numGroups x 1) cell array; C.moments{groupIdx) is a
%                       (yDims(groupIdx) x 1) cell array, and each element
%                       is a (xDim x xDim) matrix giving the posterior
%                       second moment of a row of C.
%         alpha.a    -- (numGroups x 1) array; shape parameters of alpha 
%                       posterior
%         alpha.b    -- (numGroups x xDim) array; scale parameters of 
%                       alpha posterior
%         alpha.mean -- (numGroups x xDim) array; mean precisions of
%                       loading weights (for ARD); alpha.a ./ alpha.b
%         phi.a      -- float; shape parameter of phi posterior
%         phi.b      -- (yDim x 1) array; scale parameters of phi posterior
%         phi.mean   -- (yDim x 1) array; mean precisions of observations; 
%                       alpha.a ./ alpha.b
%         xDim       -- int; number of latent variables
%         yDims      -- (1 x numGroups) array; dimensionalities of each 
%                       observed group
%     cutoff_sharedvar -- (1 x numThresh) array; list of shared variance
%                         thresholds to search over (in range [0 1]).
%     cutoff_snr       -- float; minimum signal-to-noise ratio that a
%                         group must have for ANY latents to be considered 
%                         significant
%     numBoot          -- int; Number of bootstrap samples (trials).
%
%     Optional:
%
%     verbose -- logical; set true to print status of which bootstrap
%                sample is currently being evaluated.
%
% Outputs:
%
%     R2.mean  -- (1 x numThresh) array; mean leave-group-out R^2 values
%                 across bootstrap samples
%     R2.sem   -- (1 x numThresh) array; SEM of R^2 across bootstrap
%                 samples
%     MSE.mean -- (1 x numThresh) array; leave-group-out mean-squared 
%                 errors across bootstrap samples.
%     MSE.sem  -- (1 x numThresh) array; SEM of MSE across bootstrap
%                 samples
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     17 Jan 2023 -- Initial full revision.

verbose = false;
extraOpts   = assignopts(who, varargin);

numThresh = length(cutoff_sharedvar);
N = length(seq);

% Initialize outputs
MSE_boot = nan(numBoot,numThresh);
R2_boot = nan(numBoot,numThresh);
for bootIdx = 1:numBoot
    if verbose
        fprintf('Bootstrap sample %d of %d...\n', bootIdx, numBoot); 
    end
    % Draw bootstrap sample
    bootSamples = datasample(1:N, N, 'Replace', true);
    seqBoot = seq(bootSamples);
    % Evaluate shared variance thresholds on the current bootstrap sample
    [R2_boot(bootIdx,:), MSE_boot(bootIdx,:)] ...
        = evalDimThresh(seqBoot, params, cutoff_sharedvar, cutoff_snr);
end

% Compute mean and SEM across bootstrap samples
R2.mean = mean(R2_boot,1);
R2.sem = std(R2_boot,0,1) ./ sqrt(numBoot);
MSE.mean = mean(MSE_boot,1);
MSE.sem = std(MSE_boot,0,1) ./ sqrt(numBoot);