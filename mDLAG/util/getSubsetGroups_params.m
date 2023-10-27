function params = getSubsetGroups_params(params, groups)
%
% params = getSubsetGroups_params(params, groups)
%
% Description: Get the mDLAG model parameters corresponding to the 
%              subset of groups in 'groups'.
%
% Arguments:
%
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
%
%     groups -- (1 x numKeptGroups) array; list of indices of the desired
%               group subset.
%
% Outputs:
%
%     params -- Same structure as 'params' above, but with only the
%               parameters corresponding to the groups in 'groups'.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     22 Oct 2022 -- Initial full revision.

% Initialize other relevant variables
yDims = params.yDims;
block_idxs = get_block_idxs(yDims);

params.D = params.D(groups,:);
params.C.means = params.C.means(groups);
if isfield(params.C, 'covs')
    params.C.covs = params.C.covs(groups);
end
params.C.moments = params.C.moments(groups);
params.alpha.a = params.alpha.a(groups);
params.alpha.b = params.alpha.b(groups,:);
params.alpha.mean = params.alpha.mean(groups,:);

keptGroupIdxs = [];
for groupIdx = groups
    currGroup = block_idxs{groupIdx};
    keptGroupIdxs = [keptGroupIdxs currGroup(1):currGroup(2)];
end

params.d.mean = params.d.mean(keptGroupIdxs);
params.d.cov = params.d.cov(keptGroupIdxs);
params.phi.mean = params.phi.mean(keptGroupIdxs);
params.phi.b = params.phi.b(keptGroupIdxs);

params.yDims = params.yDims(groups);