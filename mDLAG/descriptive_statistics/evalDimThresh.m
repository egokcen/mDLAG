function [R2, MSE] = evalDimThresh(seq, params, cutoff_sharedvar, cutoff_snr)
%
% [R2, MSE] = evalDimThresh(seq, params, cutoff_sharedvar, cutoff_snr)
%
% Description: Evaluate leave-group-out predictive performance across a
%              specified range of shared variance cutoffs. For each cutoff,
%              if a dimension does not exceed the treshold then it will
%              be zeroed out during prediction.
% 
% Arguments:
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
%
% Outputs:
%
%     R2      -- (1 x numThresh) array; leave-group-out R^2 values
%     MSE     -- (1 x numThresh) array; leave-group-out mean-squared errors
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     17 Jan 2023 -- Initial full revision.

numThresh = length(cutoff_sharedvar);
numGroups = length(params.yDims);

% Initialize outputs
MSE = nan(1,numThresh);
R2 = nan(1,numThresh);
for j = 1:numThresh
    % Determine significant dimensions
    [~,sigdims,~,~] = computeDimensionalities(params, ...
                                              cutoff_sharedvar(j), ...
                                              cutoff_snr);
    % Zero-out the posterior mean and second moment of C for dimensions
    % that fall below the cutoff
    sigparams = params;
    for groupIdx = 1:numGroups
        sigparams.C.means{groupIdx}(:,~sigdims(groupIdx,:)) = 0;
        for yIdx = 1:params.yDims(groupIdx)
            sigparams.C.moments{groupIdx}{yIdx}(~sigdims(groupIdx,:),:) = 0;
            sigparams.C.moments{groupIdx}{yIdx}(:,~sigdims(groupIdx,:)) = 0;
        end
    end
    % Evaluate leave-group-out predictive performance
    [R2(j), MSE(j)] = pred_mdlag(seq, sigparams);
end