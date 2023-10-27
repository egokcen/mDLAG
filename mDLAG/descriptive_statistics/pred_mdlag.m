function [R2, MSE] = pred_mdlag(seq, params)
%
% [R2, MSE] = pred_mdlag(seq, params)
%
% Description: Performs leave-group-out prediction using an existing mDLAG
%              model.
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
%
% Outputs:
%
%     MSE     -- float; leave-group-out mean-squared error
%     R2      -- float; leave-group-out R^2
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     22 Oct 2022 -- Initial full revision.

yDims = params.yDims;
xDim = params.xDim;
numGroups = length(yDims);
block_idxs = get_block_idxs(yDims);

Ys_true = seq2cell2D(seq, yDims, 'datafield', 'y');

% Construct joint data matrix
Ytrue = cat(1,Ys_true{:});
[yDim, NT] = size(Ytrue);

Ys_pred = cell(1,numGroups);
for groupIdx = 1:numGroups
    Ys_pred{groupIdx} = nan(size(Ys_true{groupIdx}));
end

% Perform leave-group-out prediction
for groupIdx = 1:numGroups
    
    targetGroup = groupIdx; % Group to be left out
    sourceGroups = setdiff(1:numGroups,targetGroup); % Observed groups
    
    % Infer latent variables given source groups
    paramsSource = getSubsetGroups_params(params,sourceGroups);
    seqSource = getSubsetGroups_seq(seq,yDims,sourceGroups);
    [seqSource,~,~] = inferX(seqSource,paramsSource);
    
    % Infer latent variables for target group, given latents for source
    % groups
    seqTarget = predX(targetGroup,params,seqSource);
    Xtarget = [seqTarget.xsm];
    
    % Predict observations for target group
    paramsTarget = getSubsetGroups_params(params,targetGroup);
    Ys_pred{targetGroup} = paramsTarget.C.means{1} * Xtarget ...
        + repmat(paramsTarget.d.mean, [1 NT]);
end

% Compute performance metrics
Ypred = cat(1,Ys_pred{:});
% MSE
MSE = immse(Ypred, Ytrue);
% R2
RSS = sum( sum( ( Ytrue - Ypred ).^2, 1 ) );
TSS = sum( sum( ( Ytrue - repmat( mean(Ytrue,2), [1 size(Ytrue,2)] ) ).^2, 1 ) );
R2 = 1 - RSS / TSS;