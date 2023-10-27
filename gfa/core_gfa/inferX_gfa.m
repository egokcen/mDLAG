function X = inferX_gfa(Ys, model)
%
% X = inferX(Ys, model)
%
% Description: Given a group factor analysis (GFA) model, infer latent
%              variables X from observations in Ys.
%
% Arguments:
%
%     Required:
%
%     Ys          -- (1 x numGroups) cell array; list of data matrices 
%                    {(y1Dim x N), (y2Dim x N), ...}
%     model       -- structure containing the following (relevant) fields:
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
%         phi.a      -- float; shape parameter of phi posterior
%         phi.b      -- (yDim x 1) array; scale parameters of phi posterior
%         phi.mean   -- (yDim x 1) array; mean precisions of observations; 
%                       alpha.a ./ alpha.b
%
% Outputs:
%     X.mean     -- (xDim x N) array; posterior mean of latent variables
%     X.cov      -- (xDim x xDim) array; posterior covariance of latent
%                   variables
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     09 Jul 2022 -- Initial full revision.

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
Call = vertcat(model.C.means{:});
[~, xDim] = size(Call);

% Zero-center observations with posterior estimate of means
Y0 = Y - repmat(model.d.mean,1,N);

% Infer latent variables, X

% Covariance
X.cov = eye(xDim);
for groupIdx = 1:numGroups
    currGroup = block_idxs{groupIdx};
    phi_m = model.phi.mean(currGroup(1):currGroup(2));
    for yIdx = 1:yDims(groupIdx)
        X.cov = X.cov + phi_m(yIdx).*model.C.moments{groupIdx}{yIdx};
    end
end
X.cov = inv(X.cov);
X.cov = 0.5 * (X.cov + X.cov'); % Ensure symmetry

% Mean
X.mean = X.cov * Call' * diag(model.phi.mean) * Y0;
    