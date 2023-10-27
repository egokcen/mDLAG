function seqPred = predX(predGroup, params, seqSource, varargin)
%
% seqPred = predX(predGroup, params, seqSource, ...)
%
% Description: Infer latent variables for the group 'predGroup' given 
%              the latent variables for the remaining groups in seqSource.
%
% Arguments:
%
%     Required:
%
%     predGroup -- int; index of group to be predicted
%     params    -- Structure containing mDLAG model parameters. 
%                  Contains the fields
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
%     seqSource -- data structure, whose nth entry (corresponding to
%                  the nth trial) has fields
%                    trialId      -- unique trial identifier
%                    T (1 x 1)    -- number of timesteps
%                    xsm (xDim*(numGroups-1) x T) -- latent variables for
%                                                    all groups but
%                                                    predGroup
%
%     Optional:
%
%     sourceGroups -- (1 x numSourceGroups) array; indices of desired
%                     groups. By default, sourceGroups is the complement
%                     of predGroup, i.e., it includes all groups BUT
%                     predGroup.
%
% Outputs:
%
%     seqPred -- data structure, whose nth entry (corresponding to
%                the nth trial) has fields
%                    trialId        -- unique trial identifier
%                    T (1 x 1)      -- number of timesteps
%                    xsm (xDim x T) -- latent variables for predGroup
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     22 Oct 2022 -- Initial full revision.

% Initialize other relevant variables
yDims = params.yDims;  
xDim = params.xDim;
numGroups = length(yDims);
sourceGroups = setdiff(1:numGroups,predGroup);
assignopts(who,varargin);

% Initialize output structure
for n = length(seqSource)
    seqPred(n).trialId = seqSource(n).trialId;
    seqPred(n).T = seqSource(n).T;
    seqPred(n).xsm = zeros(xDim,seqSource(n).T);
end

% Group trials of same length together
Tall = [seqSource.T];
Tu = unique(Tall);

% Overview:
% - Outer loop on each element of Tu.
% - For each element of Tu, find all trials with that length.
% - Do inference for all those trials together.

for j = 1:length(Tu)
    T = Tu(j);
    
    % Construct predictive GP covariance matrix
    Kpred = construct_K_pred(params,T,predGroup,'sourceGroups',sourceGroups);
    
    % Process all trials with length T
    nList    = find(Tall == T);

    % Predict latent variables for predGroup from the latent variables in 
    % sourceGroups
    for n = nList
      xsm_pred = Kpred * reshape(seqSource(n).xsm, [], 1);
      seqPred(n).xsm = reshape(xsm_pred, xDim, T);
    end
    
end
