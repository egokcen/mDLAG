function [seq, params, lbterms] = inferX(seq, params)
%
% [seq, params, lbterms] = inferX(seq, params)
%
% Description: Infer latent variables X given observations Y (in seq)
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
%
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
%     seq    -- data structure whose nth entry has new fields (these fields 
%               are added to existing fields in the seq input argument)
%                 xsm -- ((numGroups*xDim) x T) array; posterior mean 
%                        at each timepoint
%     params -- Structure containing mDLAG model parameters with new field
%                 X.cov -- data structure whose jth entry, corresponding to a group of
%                          trials of the same length, has fields
%                            T     -- int; number of time steps for this
%                                     trial group
%                            Vsm   -- (xDim*numGroups x xDim*numGroups x T)
%                                     array; posterior covariance at each 
%                                     timepoint
%                            VsmGP -- (numGroups*T x numGroups*T x xDim) 
%                                     array; posterior covariance of each 
%                                     GP
%               NOTE: For mDLAG, posterior covariances of X are the same 
%                     for trials of the same length.
%     lbterms -- Structure containing terms that contribute to the
%                variational lower bound:
%                  logdet_SigX  -- float; log|SigX|
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     13 Oct 2022 -- Initial full revision.

% Initialize other relevant variables
yDims = params.yDims;  
xDim = params.xDim;
numGroups = length(yDims);
block_idxs = get_block_idxs(yDims);

% Total number of latents, for all groups
mp = numGroups * xDim;

% Compute C'*Phi, an intermediate term
CPhi = blkdiag(params.C.means{:})' * diag(params.phi.mean);

% Compute the expectation of C'*Phi*C, an intermediate term
CPhiC = cell(1,numGroups);
for groupIdx = 1:numGroups
    currGroup = block_idxs{groupIdx};
    phi_m = params.phi.mean(currGroup(1):currGroup(2));
    CPhiC{groupIdx} = zeros(xDim,xDim);
    for yIdx = 1:yDims(groupIdx)
        CPhiC{groupIdx} = CPhiC{groupIdx} + phi_m(yIdx).*params.C.moments{groupIdx}{yIdx};
    end
end
% Collect CPhiC from each group into a large 
% (numGroups*xDim x numGroups*xDim) matrix
CPhiC = blkdiag(CPhiC{:});

% Group trials of same length together
Tall = [seq.T];
Tu = unique(Tall);

% Overview:
% - Outer loop on each element of Tu.
% - For each element of Tu, find all trials with that length.
% - Do inference for all those trials together.
lbterms.logdet_SigX = 0;
for j = 1:length(Tu)
    T = Tu(j);
    params.X.cov(j).T = T;
    
    % Construct K_big
    K_big = construct_K_mdlag(params,T);
    try
        K_big_inv = invChol_mex(K_big);
    catch
        K_big_inv = inv(K_big);
    end
    
    % Compute the full posterior covariance of X
    CPhiC_big = cell(1,T);
    [CPhiC_big{:}] = deal(CPhiC);
    try
        SigX = invChol_mex(K_big_inv + blkdiag(CPhiC_big{:}));
    catch
        SigX = inv(K_big_inv + blkdiag(CPhiC_big{:}));
    end
    
    % Now, take only subsets of SigX that will have future use
    
    % (xDim*numGroups) X (xDim*numGroups) Posterior covariance for each timepoint
    params.X.cov(j).Vsm = nan(mp,mp,T);
    idx = 1:mp;
    for t = 1:T
        cIdx = mp*(t-1)+idx;
        params.X.cov(j).Vsm(:,:,t) = SigX(cIdx,cIdx);
    end
    
    % (numGroups*T) x (numGroups*T) Posterior covariance for each GP
    
    % Posterior covariance for across-group
    params.X.cov(j).VsmGP = nan(numGroups*T,numGroups*T,xDim);
    idxs = extractLatentIndices(xDim,T,numGroups);
    for i = 1:xDim
        idx = idxs{i};
        params.X.cov(j).VsmGP(:,:,i) = SigX(idx,idx);
    end
    
    % Process all trials with length T
    nList    = find(Tall == T);
    Y0 = bsxfun(@minus, [seq(nList).y], params.d.mean); % yDim x length(nList)*T
    CPhiY0 = reshape(CPhi * Y0, mp*T, []); % (xDim*numGroups*T) x length(nList)

    % Posterior mean of X
    xsmMat  = SigX * CPhiY0; % (xDim*numGroups*T) x length(nList)
    % Reshape X and collect it in seq
    ctr = 1;
    for n = nList
      seq(n).xsm   = reshape(xsmMat(:,ctr), mp, T);      
      ctr = ctr + 1;
    end
    
    % Lower bound term
    lbterms.logdet_SigX = lbterms.logdet_SigX + length(nList)*logdet(SigX);
    
end
