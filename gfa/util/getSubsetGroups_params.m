function outparams = getSubsetGroups_params(inparams, groups)
%
% outparams = getSubsetGroups_params(inparams, groups)
%
% Description: Get the GFA model parameters corresponding to the 
%              subset of groups in 'groups'.
%
% Arguments:
%
%     inparams  -- Structure containing GFA model parameters. 
%                  Contains the fields
%         xDim       -- int; number of factors
%         yDims      -- (1 x numGroups) array; dimensionalities of each 
%                       observed group
%         X.mean     -- (xDim x N) array; posterior mean of latent
%                       variables
%         X.cov      -- (xDim x xDim) array; posterior covariance of latent
%                       variables
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
%         alpha.a    -- (numGroups x 1) array; shape parameters of 
%                       alpha posterior
%         alpha.b    -- (numGroups x xDim) array; scale parameters of 
%                       alpha posterior
%         alpha.mean -- (numGroups x xDim) array; mean precisions of
%                       loading weights (for ARD); alpha.a ./ alpha.b
%         phi.a      -- float; shape parameter of phi posterior
%         phi.b      -- (yDim x 1) array; scale parameters of phi posterior
%         phi.mean   -- (yDim x 1) array; mean precisions of observations; 
%                       alpha.a ./ alpha.b
%
%     groups -- (1 x numKeptGroups) array; list of indices of the desired
%               group subset.
%
% Outputs:
%
%     outparams -- Same structure as 'inparams' above, but with only the
%                  parameters corresponding to the groups in 'groups'.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Nov 2022 -- Initial full revision.

% Initialize other relevant variables
yDims = inparams.yDims;
block_idxs = get_block_idxs(yDims);

if isfield(inparams, 'X')
    outparams.X = inparams.X; 
end

outparams.C.means = inparams.C.means(groups);
if isfield(inparams.C, 'covs')
    outparams.C.covs = inparams.C.covs(groups);
end
outparams.C.moments = inparams.C.moments(groups);
outparams.alpha.a = inparams.alpha.a(groups);
outparams.alpha.b = inparams.alpha.b(groups,:);
outparams.alpha.mean = inparams.alpha.mean(groups,:);

keptGroupIdxs = [];
for groupIdx = groups
    currGroup = block_idxs{groupIdx};
    keptGroupIdxs = [keptGroupIdxs currGroup(1):currGroup(2)];
end

outparams.d.mean = inparams.d.mean(keptGroupIdxs);
outparams.d.cov = inparams.d.cov(keptGroupIdxs);
outparams.phi.mean = inparams.phi.mean(keptGroupIdxs);
outparams.phi.b = inparams.phi.b(keptGroupIdxs);
outparams.phi.a = inparams.phi.a;

outparams.yDims = inparams.yDims(groups);
outparams.xDim = inparams.xDim;