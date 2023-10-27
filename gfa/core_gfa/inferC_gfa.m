function model = inferC_gfa(Ys, X, model)
%
% model = inferC_gfa(Ys, model)
%
% Description: Given a group factor analysis (GFA) model, latent
%              variable estimates X, and observations in Ys, infer the
%              loading matrix C.
%
% Arguments:
%
%     Required:
%
%     Ys          -- (1 x numGroups) cell array; list of data matrices 
%                    {(y1Dim x N), (y2Dim x N), ...}
%     X.mean      -- (xDim x N) array; posterior mean of latent variables
%     X.cov       -- (xDim x xDim) array; posterior covariance of latent
%                    variables
%     model       -- structure containing the following (relevant) fields:
%         d.mean     -- (yDim x 1) array; posterior mean of mean parameter
%         d.cov      -- (yDim x 1) array; diagonal elements of the
%                       posterior covariance matrix of d
%         phi.a      -- float; shape parameter of phi posterior
%         phi.b      -- (yDim x 1) array; scale parameters of phi posterior
%         phi.mean   -- (yDim x 1) array; mean precisions of observations; 
%                       alpha.a ./ alpha.b
%
% Outputs:
%     C.means    -- (numGroups x 1) cell array; yDims(groupIdx) x xDim
%                   mean loadings matrix for each group
%     C.covs     -- (numGroups x 1) cell array; C.covs{groupIdx) is a
%                   (yDims(groupIdx) x 1) cell array, and each element
%                   is a (xDim x xDim) matrix giving the posterior
%                   covariance of a row of C.
%     C.moments  -- (numGroups x 1) cell array; C.moments{groupIdx) is a
%                   (yDims(groupIdx) x 1) cell array, and each element
%                   is a (xDim x xDim) matrix giving the posterior
%                   second moment of a row of C.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     26 Oct 2022 -- Initial full revision.

numGroups = length(Ys);
% Get individual group sizes
yDims = zeros(1,numGroups);
for groupIdx = 1:numGroups
    yDims(groupIdx) = size(Ys{groupIdx},1); 
end
block_idxs = get_block_idxs(yDims);

% Construct joint data, loadings matrices
Y = cat(1,Ys{:});
[~, N] = size(Y);
[xDim, ~] = size(X.mean);

% Zero-center observations with posterior estimate of means
Y0 = Y - repmat(model.d.mean,1,N);

% Correlation matrix between current estimate of latent variables and
% zero-centered observations.
XY = X.mean*Y0';   % (xDim x yDim) array
% Second moment
XX = N.*X.cov + X.mean*X.mean';

% Allocate output structures
model.C.means = cell(1,numGroups);
model.C.covs = cell(1,numGroups);
model.C.moments = cell(1,numGroups);
for groupIdx = 1:numGroups
    model.C.means{groupIdx} = zeros(yDims(groupIdx),xDim);
    model.C.covs{groupIdx} = cell(yDims(groupIdx),1);
    model.C.moments{groupIdx} = cell(yDims(groupIdx),1);
    for yIdx = 1:yDims(groupIdx)
        model.C.covs{groupIdx}{yIdx} = zeros(xDim);
        model.C.moments{groupIdx}{yIdx} = zeros(xDim);
    end
end

% Infer C
for groupIdx = 1:numGroups
    currGroup = block_idxs{groupIdx};
    phi_m = model.phi.mean(currGroup(1):currGroup(2));
    alpha_m = diag(model.alpha.mean(groupIdx,:));  % (xDim x xDim) array
    XY_m = XY(:,currGroup(1):currGroup(2));
    for yIdx = 1:yDims(groupIdx)
        % Covariance
        covC = inv(alpha_m + phi_m(yIdx) .* XX);
        covC = 0.5 * (covC + covC'); % Ensure symmetry
        model.C.covs{groupIdx}{yIdx} = covC;
        % Mean
        model.C.means{groupIdx}(yIdx,:) = (phi_m(yIdx) * covC * XY_m(:,yIdx))';
        % Second moment
        model.C.moments{groupIdx}{yIdx} = covC ...
            + model.C.means{groupIdx}(yIdx,:)' * model.C.means{groupIdx}(yIdx,:);
    end
end
    