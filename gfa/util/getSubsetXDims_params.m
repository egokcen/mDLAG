function outparams = getSubsetXDims_params(inparams, xDims)
%
% outparams = getSubsetXDims_params(inparams, xDims)
%
% Description: Return a GFA model with only parameters corresponding to
%              the specified latents.
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
%     xDims     -- (1 x numDims) array; latent state dimensions to be  
%                  retained in outparams.
%
% Outputs:
%
%     outparams -- Structure containing subset of mDLAG model parameters.
%                  Same format as inparams, above.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Nov 2022 -- Initial full revision.

numGroups = length(inparams.yDims);

outparams.xDim = length(xDims);
outparams.yDims = inparams.yDims;

% Latent variables, X
if isfield(inparams, 'X')
    outparams.X.mean = inparams.X.mean(xDims,:);
    outparams.X.cov = inparams.X.cov(xDims,xDims);
end

% Mean offset, d
outparams.d = inparams.d;

% ARD parameters
outparams.alpha.mean = inparams.alpha.mean(:,xDims);
outparams.alpha.b = inparams.alpha.b(:,xDims);
outparams.alpha.a = inparams.alpha.a;

% Loading matrix, C
for groupIdx = 1:numGroups
    outparams.C.means{groupIdx} = inparams.C.means{groupIdx}(:,xDims);
    for yIdx = 1:outparams.yDims(groupIdx)
        if isfield(inparams.C, 'covs')
            outparams.C.covs{groupIdx}{yIdx} ...
                = inparams.C.covs{groupIdx}{yIdx}(xDims,xDims);
        end
        outparams.C.moments{groupIdx}{yIdx} ...
            = inparams.C.moments{groupIdx}{yIdx}(xDims,xDims);
    end
end

% Observation noise, phi
outparams.phi = inparams.phi;
