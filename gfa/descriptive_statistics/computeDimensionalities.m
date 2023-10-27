function [dims, sigDims, varExp, dimTypes] = computeDimensionalities(model, cutoff_sharedvar, cutoff_snr)
%
% [dims, sigDims, varExp, dimTypes] = computeDimensionalities(model, cutoff_sharedvar, cutoff_snr)
%
% Description: Compute the number of each possible type of dimension, along
%              with the shared variance in each group explained by each
%              type of dimension.
% 
% Arguments:
%
%     model -- structure containing the following relevant GFA model
%              parameters:
%         alpha.mean -- (numGroups x xDim) array; mean precisions of
%                       loading weights
%         C.moments  -- (numGroups x 1) cell array; C.moments{groupIdx) is 
%                       a (yDims(groupIdx) x 1) cell array, and each element
%                       is a (xDim x xDim) matrix giving the posterior
%                       second moment of a row of C.
%         phi.mean   -- (yDim x 1) array; mean precisions of observations
%
%     cutoff_sharedvar -- float in [0,1]; minimum fraction of shared
%                         variance within a group that must be explained by 
%                         a latent to be considered significant.
%     cutoff_snr       -- float; minimum signal-to-noise ratio that a
%                         group must have for ANY latents to be considered 
%                         significant
%
% Outputs:
%
%     dims     -- (1 x numDimTypes) array; The number of each type of
%                 dimension. dims(i) corresponds to the dimension type in 
%                 dimTypes(:,i).
%     sigDims  -- (numGroups x xDim) logical array; sigDims(i,j) is 1 if
%                 latent j explains significant shared variance in group i,
%                 0 otherwise.
%     varExp   -- (numGroups x numDimTypes) array; varExp(i,j) is the
%                 proportion of shared variance explained in group i by
%                 dimension type j. varExp(:,j) corresponds to the
%                 dimension type in dimTypes(:,j).
%     dimTypes -- (numGroups x numDimTypes) array; dimTypes(:,i) is a
%                 binary vector indicating the structure of dimension type
%                 i. 1 indicates a group is involved, 0 otherwise.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     07 Sep 2022 -- Initial full revision.

numGroups = size(model.alpha.mean,1);

% Determine all dimension types
dimTypes = generateDimTypes(numGroups);
numDimTypes = size(dimTypes,2);

% Compute signal-to-noise ratios
snr = computeSNR(model.C, model.phi);

% Relative shared variance explained by each dimension
alpha_inv = 1./model.alpha.mean;
alpha_inv_rel = alpha_inv ./ sum(alpha_inv,2);

% Take dimensions only if shared variance and SNR are significant
sigDims = alpha_inv_rel > cutoff_sharedvar & snr > cutoff_snr;
dims = nan(1,numDimTypes);
varExp = nan(numGroups,numDimTypes);
for dimIdx = 1:numDimTypes 
    dims(dimIdx) = sum(ismember(sigDims',dimTypes(:,dimIdx)','rows'));
    varExp(:,dimIdx) = sum(alpha_inv_rel(:,ismember(sigDims',dimTypes(:,dimIdx)','rows')),2);
end
