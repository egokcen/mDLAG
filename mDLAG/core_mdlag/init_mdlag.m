function params = init_mdlag(seq,yDims,xDim,binWidth,varargin)
%
% params = init_mdlag(seq,yDims,xDim,binWidth,...)
%
% Description: Initialize mDLAG parameters for the variational EM algorithm.
%
% Arguments:
%
%     Required:
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- observed data
%     yDims    -- (1 x numGroups) array; dimensionalities of each observed 
%                 group
%     xDim     -- int; number of latent variables
%     binWidth -- float; bin width or sample period, in units of time.
%                 Note: For all other temporal variables, keep units
%                       consistent with binWidth.
%
%     Optional:
%    
%     startTau   -- float; Initial GP timescale, in units oftime 
%                   (default: 2*binWidth)
%     startEps   -- float; Initial GP noise variance (default: 1e-3)
%     startDelay -- (numGroups x xDim) array; Initial delays between 
%                   latents and groups. Entries in units of time. 
%                   (default: [])
%     covType    -- string; Specify GP covariance kernel type. Options
%                   currently supported:
%                       'rbf' -- Radial basis function, or squared
%                                exponential kernel
%                   (default: 'rbf')
%     prior      -- structure with the following fields:
%                    d.beta  -- positive float; precision of mean parameter
%                               generative model (Gaussian)
%                    phi.a   -- positive float; 'a' shape parameter of
%                               observation precision (phi) generative 
%                               model (Gamma with mean a/b)
%                    phi.b   -- positive float; 'b' scale parameter of
%                               observation precision (phi) generative
%                               model (Gamma with mean a/b)
%                    alpha.a -- positive float; 'a' shape parameter of 
%                               alpha parameter generative model 
%                               (Gamma with mean a/b)
%                    alpha.b -- positive float; 'b' scale parameter of
%                               alpha parameter generative model 
%                               (Gamma with mean a/b)
%                    (default: 1e-12 for all values)
%     randomSeed -- int or string; seed the random number generator, 
%                   for reproducible initializations (default: 'shuffle')
%     verbose    -- boolean; specifies whether to display status messages
%                   (default: false)
%     saveCcov   -- logical; set true to save posterior covariance and 
%                   of C. For large yDim and xDim, these structures can use
%                   a lot of memory. (default: false)
% 
% Outputs:
%
%     params -- Structure containing mDLAG model parameters at which EM
%               algorithm can be initialized. Contains the fields
% 
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
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     05 Oct 2022 -- Initial full revision.
%     23 Oct 2022 -- Added the performance-enhancing saveCcov option.

startTau      = 2*binWidth; % in units of time
startDelay    = [];         % in units of time
startEps      = 1e-3;
covType       = 'rbf';
prior.d.beta = 1e-12;
prior.phi.a = 1e-12;
prior.phi.b = 1e-12;
prior.alpha.a = 1e-12;
prior.alpha.b = 1e-12;
randomSeed = 'shuffle';
verbose = false;
saveCcov = false;
extraOpts     = assignopts(who, varargin);

% Seed random number generator
rng(randomSeed);

numGroups = length(yDims);
params.xDim = xDim;
params.yDims = yDims;
block_idxs = get_block_idxs(yDims);

% Concatenate data across trials
Ys = seq2cell2D(seq, yDims, 'datafield', 'y');

% Construct joint data matrix
Y = cat(1,Ys{:});
[yDim, NT] = size(Y);
covY = cov(Y', 1);

% =======================================
% Initialize observation model parameters
% =======================================

% Mean parameter
params.d.mean = mean(Y,2);                      % Mean
params.d.cov = (1/prior.d.beta).*ones(yDim,1);  % Covariance

% Noise precisions
params.phi.a = prior.phi.a + NT/2;
params.phi.b = prior.phi.b .* ones(yDim,1);
params.phi.mean = 1./diag(covY);

% Loading matrices
params.C.means = cell(1,numGroups);
if saveCcov
    params.C.covs = cell(1,numGroups);
end
params.C.moments = cell(1,numGroups);
for groupIdx = 1:numGroups
    params.C.means{groupIdx} = zeros(yDims(groupIdx),xDim);
    covY_m = cov(Ys{groupIdx}', 1);
    if rank(covY_m) == yDims(groupIdx)
        scale = exp(2*sum(log(diag(chol(covY_m))))/yDim); % Scale by determinant of covY_m
    else
        % covY_m may not be full rank because N < yDim
        if verbose
            fprintf('WARNING in init_mdlag: Data matrix for group %d is not full rank.\n', groupIdx);
        end
        r     = rank(covY_m);
        e     = sort(eig(covY_m), 'descend');
        scale = geomean(e(1:r));
    end
    params.C.means{groupIdx} = randn(yDims(groupIdx),xDim)*sqrt(scale/xDim);
    if saveCcov
        params.C.covs{groupIdx} = cell(yDims(groupIdx),1);
    end
    params.C.moments{groupIdx} = cell(yDims(groupIdx),1);
    for yIdx = 1:yDims(groupIdx)
        covC = zeros(xDim);
        if saveCcov
            params.C.covs{groupIdx}{yIdx} = covC;
        end
        params.C.moments{groupIdx}{yIdx} = covC ...
            + params.C.means{groupIdx}(yIdx,:)'*params.C.means{groupIdx}(yIdx,:);
    end
end

% ARD parameters
params.alpha.a = prior.alpha.a + (yDims./2)';  % (numGroups x 1) array
params.alpha.b = prior.alpha.b .* ones(numGroups,xDim);
params.alpha.mean = nan(numGroups,xDim);
% Scale ARD parameters to match the data
for groupIdx = 1:numGroups
    CC_m = zeros(xDim);
    for yIdx = 1:yDims(groupIdx)
        CC_m = CC_m + params.C.moments{groupIdx}{yIdx}; 
    end
    params.alpha.mean(groupIdx,:) = (yDims(groupIdx)./diag(CC_m))';
end

% ======================================
% Initialize Gaussian process parameters
% ======================================

params.covType = covType;

% Delay matrix
if isempty(startDelay)
    % No initial delays were specified. Initialize to zeros.
    params.D = zeros(numGroups,xDim);
else
    % Assign initial delays to the specified values
    assert(isequal(size(startDelay),[numGroups xDim]), ... 
           'startDelay must have size [numGroups xDim]');
    params.D = startDelay;
end

% GP timescale
% Assume binWidth is the time step size.
params.gamma = (binWidth ./ startTau).^2 .* ones(1, xDim);
% GP noise variance
params.eps = startEps .* ones(1, xDim);