function params = inferC(seq, params)
%
% params = inferC(seq, params)
%
% Description: Infer loading matrix C given observations Y (in seq)
%              and mDLAG model parameters in params.
%
% Arguments:
%
%     Required:
%
%     seq     -- data structure, whose nth entry (corresponding to
%                the nth trial) has fields
%                    trialId      -- unique trial identifier
%                    T (1 x 1)    -- number of timesteps
%                    y (yDim x T) -- neural data
%                    xsm          -- ((numGroups*xDim) x T) array; 
%                                    posterior mean at each timepoint
%
%     params  -- Structure containing mDLAG model parameters. 
%                Contains the fields
%         covType    -- string; type of GP covariance (e.g., 'rbf')
%         X.cov      -- data structure whose jth entry, corresponding
%                       to a group of trials of the same length, has fields
%                         T     -- int; number of time steps for this
%                                  trial group
%                         Vsm   -- (xDim*numGroups x xDim*numGroups x T)
%                                  array; posterior covariance at each 
%                                  timepoint
%                         VsmGP -- (numGroups*T x numGroups*T x xDim) 
%                                  array; posterior covariance of each GP
%         gamma      -- (1 x xDim) array; GP timescales in units of 
%                       time are given by 'binWidth ./ sqrt(gamma)'                                                    
%         eps        -- (1 x xDim) array; GP noise variances
%         D          -- (numGroups x xDim) array; delays from latents to 
%                       observed variables. NOTE: Delays are reported as 
%                       (real-valued) number of time-steps.
%         d.mean     -- (yDim x 1) array; posterior mean of mean parameter
%         d.cov      -- (yDim x 1) array; diagonal elements of the
%                       posterior covariance matrix of d
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
%     params -- Structure containing mDLAG model parameters with new fields
%         C.means   -- (numGroups x 1) cell array; yDims(groupIdx) x xDim
%                      mean loadings matrix for each group
%         C.covs    -- (numGroups x 1) cell array; C.covs{groupIdx) is a
%                      (yDims(groupIdx) x 1) cell array, and each element
%                      is a (xDim x xDim) matrix giving the posterior
%                      covariance of a row of C.
%         C.moments -- (numGroups x 1) cell array; C.moments{groupIdx) is a
%                      (yDims(groupIdx) x 1) cell array, and each element
%                      is a (xDim x xDim) matrix giving the posterior
%                      second moment of a row of C.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Oct 2022 -- Initial full revision.

% Initialize relevant variables
yDims = params.yDims;
numGroups = length(yDims);
block_idxs = get_block_idxs(yDims);

% Zero-centered observations
Y0 = bsxfun(@minus, [seq.y], params.d.mean);
% Posterior means of X for each group
Xs = seq2cell2D(seq, params.xDim.*ones(1,numGroups), 'datafield', 'xsm');
% Posterior second moments of X for each group
XX = computeXmoment(seq,params);

% Allocate output structures
params.C.means = cell(1,numGroups);
params.C.covs = cell(1,numGroups);
params.C.moments = cell(1,numGroups);
for groupIdx = 1:numGroups
    params.C.means{groupIdx} = zeros(yDims(groupIdx),params.xDim);
    params.C.covs{groupIdx} = cell(yDims(groupIdx),1);
    params.C.moments{groupIdx} = cell(yDims(groupIdx),1);
    for yIdx = 1:yDims(groupIdx)
        params.C.covs{groupIdx}{yIdx} = zeros(params.xDim);
        params.C.moments{groupIdx}{yIdx} = zeros(params.xDim);
    end
end

% Infer C
for groupIdx = 1:numGroups
    currGroup = block_idxs{groupIdx};
    phi_m = params.phi.mean(currGroup(1):currGroup(2));
    alpha_m = diag(params.alpha.mean(groupIdx,:));  % (xDim x xDim) array
    XY_m = Xs{groupIdx} * Y0(currGroup(1):currGroup(2),:)';
    for yIdx = 1:yDims(groupIdx)
        % Covariance
        params.C.covs{groupIdx}{yIdx} = inv(alpha_m + phi_m(yIdx) .* XX{groupIdx});
        params.C.covs{groupIdx}{yIdx} ...
            = 0.5 * (params.C.covs{groupIdx}{yIdx} + params.C.covs{groupIdx}{yIdx}'); % Ensure symmetry
        % Mean
        params.C.means{groupIdx}(yIdx,:) ...
            = (phi_m(yIdx) * params.C.covs{groupIdx}{yIdx} * XY_m(:,yIdx))';
        % Second moment
        params.C.moments{groupIdx}{yIdx} ...
            = params.C.covs{groupIdx}{yIdx} ...
            + params.C.means{groupIdx}(yIdx,:)' * params.C.means{groupIdx}(yIdx,:);
    end
end
