function [R2, MSE] = bootstrapDimThresh(Ys, params, cutoff_sharedvar, cutoff_snr, numBoot, varargin)
%
% [R2, MSE] = bootstrapDimThresh(Ys, params, cutoff_sharedvar, cutoff_snr, numBoot, ...)
%
% Description: Evaluate leave-group-out predictive performance across a
%              specified range of shared variance cutoffs using a
%              non-parametric bootstrap procedure.
% 
% Arguments:
%
%     Required:
%
%     Ys         -- (1 x numGroups) cell array; list of data matrices 
%                   {(y1Dim x N), (y2Dim x N), ...}
%     params     -- structure containing GFA model parameters
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
%     18 Jan 2023 -- Initial full revision.

verbose = false;
extraOpts   = assignopts(who, varargin);

numThresh = length(cutoff_sharedvar);
numGroups = length(Ys);
N = size(Ys{1},2);

% Initialize outputs
MSE_boot = nan(numBoot,numThresh);
R2_boot = nan(numBoot,numThresh);
for bootIdx = 1:numBoot
    if verbose
        fprintf('Bootstrap sample %d of %d...\n', bootIdx, numBoot); 
    end
    % Draw bootstrap sample
    bootSamples = datasample(1:N, N, 'Replace', true);
    Ys_boot = cell(1,numGroups);
    for groupIdx = 1:numGroups
        Ys_boot{groupIdx} = Ys{groupIdx}(:,bootSamples); 
    end
    % Evaluate shared variance thresholds on the current bootstrap sample
    [R2_boot(bootIdx,:), MSE_boot(bootIdx,:)] ...
        = evalDimThresh(Ys_boot, params, cutoff_sharedvar, cutoff_snr);
end

% Compute mean and SEM across bootstrap samples
R2.mean = mean(R2_boot,1);
R2.sem = std(R2_boot,0,1) ./ sqrt(numBoot);
MSE.mean = mean(MSE_boot,1);
MSE.sem = std(MSE_boot,0,1) ./ sqrt(numBoot);