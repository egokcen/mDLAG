function [R2, MSE] = pred_gfa_pairwise(Ys, C, d, phi)
%
% [R2, MSE] = pred_gfa_pairwise(Ys, C, d, phi)
%
% Description: Performs prediction between each pair of groups using an 
%              existing GFA model.
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
%     R2.uni  -- (numGroups x numGroups) array; R2(i,j) gives the R^2 
%                 value when predicting group j from group i
%                 (unidirectional prediction, asymmetric)
%     R2.agg  -- (numGroups x numGroups) array; R2(i,j) gives the aggregate
%                R^2 value when predicting between groups i and j
%                (symmetric; only upper triangular portion is filled in)
%     MSE.uni -- (numGroups x numGroups) array; MSE(i,j) gives the
%                mean-squared error when predicting group j from group i
%                (unidirectional prediction, asymmetric)
%     MSE.agg -- (numGroups x numGroups) array; MSE(i,j) gives the aggregate
%                MSE value when predicting between groups i and j
%                (symmetric; only upper triangular portion is filled in)
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     08 Nov 2022 -- Initial full revision.
%     16 Nov 2022 -- Overhauled to include aggregate metrics.

numGroups = length(Ys);
pairs = nchoosek(1:numGroups,2);
numPairs = size(pairs,1);

N = size(Ys{1},2);
yDims = nan(1,numGroups);
for groupIdx = 1:numGroups
    yDims(groupIdx) = size(Ys{groupIdx},1);
end
block_idxs = get_block_idxs(yDims);
xDim = size(C.means{1},2);

% Perform pairwise prediction
R2.uni = nan(numGroups,numGroups);
R2.agg = nan(numGroups,numGroups);
MSE.agg = nan(numGroups,numGroups);
MSE.uni = nan(numGroups,numGroups);

for pairIdx = 1:numPairs
    
    pair = pairs(pairIdx,:);
    Ys_pred = cell(1,length(pair));
    for targetGroup = pair
        targetDims = block_idxs{targetGroup};
        sourceGroup = setdiff(pair, targetGroup); % Observed groups
    
        % Infer latent variables given source group
        sourceDims = block_idxs{sourceGroup};
        
        % Covariance
        X.cov = eye(xDim);
        phi_m = phi.mean(sourceDims(1):sourceDims(2));
        for yIdx = 1:yDims(sourceGroup)
            X.cov = X.cov + phi_m(yIdx).*C.moments{sourceGroup}{yIdx};
        end
        X.cov = inv(X.cov);
        X.cov = 0.5 * (X.cov + X.cov'); % Ensure symmetry

        % Mean
        X.mean = zeros(xDim,N);
        phi_m = phi.mean(sourceDims(1):sourceDims(2));
        d_m = d.mean(sourceDims(1):sourceDims(2));
        X.mean = X.mean + C.means{sourceGroup}' * diag(phi_m) ...
            * (Ys{sourceGroup} - repmat(d_m, [1 N]));
        X.mean = X.cov * X.mean;

        % Predict observations for target group
        Ypred = C.means{targetGroup} * X.mean ...
            + repmat(d.mean(targetDims(1):targetDims(2)), [1 N]);
        Ys_pred{targetGroup} = Ypred;

        % Compute performance metrics for this target group
        % MSE
        MSE.uni(sourceGroup,targetGroup) = immse(Ypred, Ys{targetGroup});
        % R2
        RSS = sum( sum( ( Ys{targetGroup} - Ypred ).^2, 1 ) );
        TSS = sum( sum( ( Ys{targetGroup} - repmat( mean(Ys{targetGroup},2), [1 size(Ys{targetGroup},2)] ) ).^2, 1 ) );
        R2.uni(sourceGroup,targetGroup) = 1 - RSS / TSS;
        
    end
    
    % Compute aggregate R^2
    Ypred = vertcat(Ys_pred{:});
    Ytrue = vertcat(Ys{pair});
    % MSE
    MSE.agg(pair(1),pair(2)) = immse(Ypred, Ytrue);
    % R2
    RSS = sum( sum( ( Ytrue - Ypred ).^2, 1 ) );
    TSS = sum( sum( ( Ytrue - repmat( mean(Ytrue,2), [1 size(Ytrue,2)] ) ).^2, 1 ) );
    R2.agg(pair(1),pair(2)) = 1 - RSS / TSS;
end