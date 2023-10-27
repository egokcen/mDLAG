function [Ds,gams,alphas] = getSubsetXDims_trackedParams(Ds,gams,alphas,xDims)
%
% [Ds,gams,alphas] = getSubsetXDims_trackedParams(Ds,gams,alphas,xDims)
%
% Description: Remove unwanted latents from mDLAG model parameters that
%              are tracked throughout fitting.
%
% Arguments:
%
%     Ds        -- (1 x numIters) cell array; the estimated delay matrix
%                  (D) after each EM iteration.
%     gams      -- (1 x numIters) cell arry; estimated gamma after each EM
%                  iteration.
%     alphas    -- (1 x numIters) cell arry; estimated ARD parameters
%                  (alpha.mean) after each EM iteration.
%     xDims     -- (1 x numDims) array; latent state dimensions to be  
%                  retained in outparams.
%
% Outputs:
%
%     Ds        -- (1 x numIters) cell array; same as Ds above, but with
%                  relevant xDims removed.
%     gams      -- (1 x numIters) cell arry; same as gams above, but with
%                  relevant xDims removed.
%     alphas    -- (1 x numIters) cell arry; same as alphas above, but with
%                  relevant xDims removed.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Oct 2022 -- Initial full revision.

numIters = length(Ds);
for i = 1:numIters
    Ds{i} = Ds{i}(:,xDims);
    gams{i} = gams{i}(xDims);
    alphas{i} = alphas{i}(:,xDims);
end
