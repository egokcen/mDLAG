function [R2, MSE] = pred_gfa(Ys, C, d, phi)
%
% [R2, MSE] = pred_gfa(Ys, C, d, phi)
%
% Description: Performs leave-group-out prediction using an existing GFA
%              model.
%
% Arguments:
%
%     Ys         -- (1 x numGroups) cell array; list of data matrices 
%                   {(y1Dim x N), (y2Dim x N), ...}
%     C.means    -- (numGroups x 1) cell array; yDims(groupIdx) x xDim
%                   mean loadings matrix for each group
%     C.moments  -- (numGroups x 1) cell array; C.moments{groupIdx) is a
%                   (yDims(groupIdx) x 1) cell array, and each element
%                   is a (xDim x xDim) matrix giving the posterior
%                   second moment of a row of C.
%     d.mean     -- (yDim x 1) array; posterior mean of mean parameter
%     phi.mean   -- (yDim x 1) array; mean precisions of observations
%
% Outputs:
%
%     R2.agg    -- float; aggregate leave-group-out R^2
%     R2.indiv  -- (1 x numGroups) array; R2.indiv(i) gives the R^2 value
%                  when predicting group i given the remaining groups
%     MSE.agg   -- float; aggregate leave-group-out mean-squared error
%     MSE.indiv -- (1 x numGroups) array; MSE.indiv(i) gives the MSE
%                  when predicting group i given the remaining groups
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     02 Jul 2022 -- Initial full revision.
%     08 Nov 2022 -- Added outputs of predictive performance for each
%                    individual target group.

numGroups = length(Ys);
N = size(Ys{1},2);
yDims = nan(1,numGroups);
Ys_pred = cell(1,numGroups);
for groupIdx = 1:numGroups
    yDims(groupIdx) = size(Ys{groupIdx},1);
    Ys_pred{groupIdx} = nan(size(Ys{groupIdx}));
end
yDim = sum(yDims);
block_idxs = get_block_idxs(yDims);
xDim = size(C.means{1},2);

% Perform leave-group-out prediction
R2.indiv = nan(1,numGroups);
MSE.indiv = nan(1,numGroups);
for groupIdx = 1:numGroups
    
    targetGroup = groupIdx; % Group to be left out
    sourceGroups = setdiff(1:numGroups, targetGroup); % Observed groups
    
    % Infer latent variables given source groups
    
    % Covariance
    X.cov = eye(xDim);
    for sourceGroup = sourceGroups
        sourceDims = block_idxs{sourceGroup};
        phi_m = phi.mean(sourceDims(1):sourceDims(2));
        for yIdx = 1:yDims(sourceGroup)
            X.cov = X.cov + phi_m(yIdx).*C.moments{sourceGroup}{yIdx};
        end
    end
    X.cov = inv(X.cov);
    X.cov = 0.5 * (X.cov + X.cov'); % Ensure symmetry
    
    % Mean
    X.mean = zeros(xDim,N);
    for sourceGroup = sourceGroups
        sourceDims = block_idxs{sourceGroup};
        phi_m = phi.mean(sourceDims(1):sourceDims(2));
        d_m = d.mean(sourceDims(1):sourceDims(2));
        X.mean = X.mean + C.means{sourceGroup}' * diag(phi_m) ...
            * (Ys{sourceGroup} - repmat(d_m, [1 N]));
    end
    X.mean = X.cov * X.mean;
    
    % Predict observations for left-out group
    targetDims = block_idxs{targetGroup};
    Ys_pred{targetGroup} = C.means{targetGroup} * X.mean ...
        + repmat(d.mean(targetDims(1):targetDims(2)), [1 N]);
    
    % Compute performance metrics for this target group
    % MSE
    MSE.indiv(targetGroup) = immse(Ys_pred{targetGroup}, Ys{targetGroup});
    % R2
    RSS = sum( sum( ( Ys{targetGroup} - Ys_pred{targetGroup} ).^2, 1 ) );
    TSS = sum( sum( ( Ys{targetGroup} - repmat( mean(Ys{targetGroup},2), [1 size(Ys{targetGroup},2)] ) ).^2, 1 ) );
    R2.indiv(targetGroup) = 1 - RSS / TSS;
end

% Compute aggregate performance metrics
Ytrue = vertcat(Ys{:});
Ypred = vertcat(Ys_pred{:});
% MSE
MSE.agg = immse(Ypred, Ytrue);
% R2
RSS = sum( sum( ( Ytrue - Ypred ).^2, 1 ) );
TSS = sum( sum( ( Ytrue - repmat( mean(Ytrue,2), [1 size(Ytrue,2)] ) ).^2, 1 ) );
R2.agg = 1 - RSS / TSS;