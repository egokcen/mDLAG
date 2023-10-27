function [currentParams,seq,trackedParams,flags] ...
    = em_mdlag(currentParams,seq,xDim,varargin)
%
% [currentParams,seq,trackedParams,flags] = em_mdlag(currentParams,seq,xDim, ...)
%
% Description: Fit a (multi-group) Delayed Latents Across Groups (mDLAG) 
%              model using a variational EM algorithm with mean-field
%              approximation.
%
% Arguments:
%
%     Required:
%
%     currentParams -- Structure containing mDLAG model parameters at which 
%                      EM algorithm is initialized. Contains the fields
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
%     seq     -- data structure, whose nth entry (corresponding to the nth 
%                trial) has fields
%                    trialId      -- unique trial identifier
%                    T (1 x 1)    -- number of timesteps
%                    y (yDim x T) -- observed data
%     xDim    -- int; number of latent variables
%
%     Optional:
%
%     prior     -- structure with the following fields:
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
%
%     tol       -- float; stopping criterion for EM (default: 1e-8)
%     maxIters  -- int; maximum number of EM iterations (default: 1e6)
%     freqLB    -- int; lower bound is computed every freqLB EM iterations.
%                  freqLB = 1 means that the lower bound is computed every
%                  iteration. (default: 10)
%     freqParam -- int; store intermediate values for delays and timescales
%                  every freqParam EM iterations (default: 100)
%     learnDelays -- logical; set true to learn delay parameters;
%                    otherwise, delays will remain fixed at their initial
%                    value.
%     verbose   -- boolean; specifies whether to display status messages
%                  (default: false)
%     minVarFrac -- float; fraction of overall data variance for each 
%                   observed dimension to set as the private variance 
%                   floor. (default: 1e-3)
%     maxDelayFrac -- float in range [0,1]; Constrain estimated delays to
%                     be no more than a certain fraction of the trial 
%                     length. (default: 0.5)
%     maxTauFrac -- float in range [0,1]; Constrain estimated timescales to
%                   be no more than a certain fraction of the trial 
%                   length. (default: 1.0)
%     trackedParams -- structure containing the tracked parameters from a
%                      previous fitting attempt, with the intention of 
%                      starting where that attempt left off. See 
%                      trackedParams below. (default: {})
%     pruneX    -- logical; set true to remove dimensions from X that
%                  become inactive in all groups. Can speed up EM runtime
%                  and improve memory efficiency. (default: true)
%     saveXcov  -- logical; set true to save posterior covariance of latent
%                  variables X. For large datasets, these structures can
%                  use a lot of memory. (default: false)
%     saveCcov  -- logical; set true to save posterior covariance and 
%                  of C. For large yDim and xDim, these structures can use
%                  a lot of memory. (default: false)
%
% Outputs:
%
%     currentParams -- Structure containing mDLAG model parameters returned 
%                  by EM algorithm (same format as currentParams) with the
%                  addition of (if saveX is true)
%                    X.cov -- data structure whose jth entry, corresponding
%                             to a group of trials of the same length, has 
%                             fields
%                               T     -- int; number of time steps for this
%                                        trial group
%                               Vsm   -- (xDim*numGroups x xDim*numGroups x T)
%                                        array; posterior covariance at each 
%                                        timepoint
%                               VsmGP -- (numGroups*T x numGroups*T x xDim) 
%                                        array; posterior covariance of each 
%                                        GP
%                             NOTE: For mDLAG, posterior covariances of X 
%                                   are the same for trials of the same 
%                                   length.
%     seq       -- data structure whose nth entry has new fields 
%                  (these fields are added to existing fields in the seq 
%                  input argument)
%                  xsm   -- ((numGroups*xDim) x T) array; posterior mean 
%                           at each timepoint
%     trackedParams -- structure containing parameters tracked throughout 
%                      fitting:
%         lb       -- (1 x numIters) array; variational lower bound at each 
%                     iteration
%         iterTime -- (1 x numIters) array; computation time for each EM
%                     iteration
%         Ds       -- (1 x numIters) cell array; the estimated delay matrix
%                     (D) after each EM iteration.
%         gams     -- (1 x numIters) cell arry; estimated gamma after each 
%                     EM iteration.
%         alphas   -- (1 x numIters) cell arry; estimated ARD parameters
%                     (alpha.mean) after each EM iteration.
%     flags     -- structure containing various status messages about
%                  the fitting procedure:
%         Convergence          -- logical; 1 if lower bound converged 
%                                 before reaching maxIters EM iterations; 
%                                 0 otherwise
%         DecreasingLowerBound -- logical; 1 if lower bound decreased 
%                                 during fitting; 0 otherwise.
%         PrivateVarianceFloor -- logical; 1 if private variance floor was 
%                                 used on any values of phi; 0 otherwise.
%         xDimsRemoved         -- int; Number of latent dimensions
%                                 removed (if pruneX is true) due to 
%                                 low variance in all groups.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     05 Oct 2022 -- Initial full revision.
%     23 Oct 2022 -- Added the performance-enhancing pruneX, saveXcov, and
%                    saveCcov options.

prior.d.beta = 1e-12;
prior.phi.a = 1e-12;
prior.phi.b = 1e-12;
prior.alpha.a = 1e-12;
prior.alpha.b = 1e-12;
tol = 1e-8; 
maxIters = 1e6;
freqLB = 10;
freqParam = 100;
learnDelays = true;
verbose = false;
minVarFrac = 1e-3;
maxDelayFrac  = 0.5;
maxTauFrac = 1.0;
trackedParams = {};
pruneX = true;
saveXcov = false;
saveCcov = false;
assignopts(who, varargin);

% Initialize common variables
yDims = currentParams.yDims;
yDim = sum(yDims);
numGroups = length(yDims);
% Useful for extracting correct-sized blocks from matrices later
block_idxs = get_block_idxs(yDims);
xDims = xDim.*ones(1,numGroups);

% If currentParams was initialized to include C.covs, but C.covs should
% not be saved, then remove it.
if ~saveCcov && isfield(currentParams.C, 'covs')
    currentParams.C = rmfield(currentParams.C, 'covs');
end

% If C.covs should be saved, but currentParams was not initialized with it,
% then initialize the structure (initial value doesn't matter here).
if saveCcov && ~isfield(currentParams.C, 'covs')
    currentParams.C.covs = currentParams.C.moments; 
end

% Concatenate data across trials
Ys = seq2cell2D(seq, yDims, 'datafield', 'y');

% Construct joint data matrix
Y = cat(1,Ys{:});
[yDim, NT] = size(Y);
Y2 = Y.^2;

% Get private variance floors
covY = cov(Y', 1);
varFloor = minVarFrac * diag(covY);

% Define constants
% Latent magnitude must be greater than this value in at least one group to 
% remain in the model
PRUNEX_TOL = 1e-7;

% Constant factors in elbo
const_lik = -((yDim*NT)/2) * log(2*pi);

% Constants in observation noise model
alogb_phi = prior.phi.a * log(prior.phi.b);
loggammaa_phi_prior = gammaln(prior.phi.a);
loggammaa_phi_post = gammaln(currentParams.phi.a);
digammaa_phi = psi(currentParams.phi.a);

% Constants in ARD model
alogb_alpha = prior.alpha.a * log(prior.alpha.b);
loggammaa_alpha_prior = gammaln(prior.alpha.a);
loggammaa_alpha_post = gammaln(currentParams.alpha.a); % (numGroups x 1) array
digammaa_alpha = psi(currentParams.alpha.a);           % (numGroups x 1) array

% Make sure initial delays are within specified constraints
% Convert maxDelayFrac to units of "time steps", the same units as in D
constraints.maxDelay = maxDelayFrac*min([seq.T]); 

% If delays are outside the range (minDelay,maxDelay), then clamp them
% within the range
currentParams.D(currentParams.D >= constraints.maxDelay) = 0.95*constraints.maxDelay;
currentParams.D(currentParams.D <= -constraints.maxDelay) = -0.95*constraints.maxDelay;

% Convert maxTauFrac to unitless quantity 'gamma'
constraints.minGamma = 1/(maxTauFrac*min([seq.T]))^2;

if isempty(trackedParams)
    % Start convergence tracking fresh
    lbi = -Inf;                   % Initial variational lower bound
    lb = [];                      % Lower bound at each iteration
    iterTime = [];                % Time it takes to complete each iteration
    Ds = {currentParams.D};       % Delays each iteration
    gams = {currentParams.gamma}; % Timescales each iteration
    alphas = {currentParams.alpha.mean}; % ARD parameters each iteration
    startIter = 1;                % Initial value for EM loop index
else
    % Start convergence tracking based on where a previous attempt left
    % off.
    lb        = trackedParams.lb;       % Lower bound at each iteration
    lbi       = trackedParams.lb(end);  % Most recent lower bound
    lbbase    = trackedParams.lb(2);    % Base lower bound
    iterTime  = trackedParams.iterTime; % Time it takes to complete each iteration
    Ds        = trackedParams.Ds;       % Delays each iteration
    gams      = trackedParams.gams;     % Timescales each iteration
    alphas    = trackedParams.alphas;   % ARD parameters each iteration
    startIter = length(lb) + 1;         % Initial value for EM loop index
end

% Initialize status flags
flags.Convergence = 0;
flags.DecreasingLowerBound = 0;
flags.PrivateVarianceFloor = 0;
flags.xDimsRemoved = 0;

% Initialize posterior mean of latent states X
[seq,~,~] = inferX(seq, currentParams);
Xs = seq2cell2D(seq, xDims, 'datafield', 'xsm');

% Begin EM iterations
for i = startIter:maxIters
    
    % Determine when to actually compute the lower bound
    if (rem(i, freqLB) == 0) || (i<=2) || (i == maxIters)
        getLB = true;
    else
        getLB = false;
    end
    
    tic; % For tracking the computation time of each iteration
    
    % ======================
    % Latent states, X
    % ======================
    
    % Check if any latent states need to be removed
    if pruneX
        % To be kept, latent states must be active in at least one group
        keptXdims = zeros(xDim,numGroups);
        for groupIdx = 1:numGroups
            keptXdims(:,groupIdx) = (mean(Xs{groupIdx}.^2,2) > PRUNEX_TOL); 
        end
        keptXdims = find(any(keptXdims,2));
        
        if length(keptXdims) < xDim
            currentParams = getSubsetXDims_params(currentParams,keptXdims);
            [Ds,gams,alphas] = getSubsetXDims_trackedParams(Ds,gams,alphas,keptXdims);
            flags.xDimsRemoved = flags.xDimsRemoved + (xDim - currentParams.xDim);
            xDim = currentParams.xDim;
            xDims = xDim.*ones(1,numGroups);
            if xDim <= 0
                % Stop fitting if no significant latent dimensions remain.
                break;
            end
        end
    end
    
    % Posterior mean and portions of the posterior covariance of X
    [seq, currentParams, lbterms] = inferX(seq, currentParams);
    
    % Posterior second moments of X for each group
    XX = computeXmoment(seq,currentParams);
    
    % Concatenate latent states across trials
    Xs = seq2cell2D(seq, xDims, 'datafield', 'xsm');
    X = cat(1,Xs{:});
    
    % ======================
    % GP parameters
    % ======================
    
    res = learnGPparams_mdlag(seq, currentParams, constraints);
    switch currentParams.covType
        case 'rbf'
            currentParams.gamma = res.gamma; 
            if learnDelays && numGroups > 1
                currentParams.D = res.D;
            end
            if (rem(i, freqParam) == 0) || (i == maxIters)
                % Store current delays and timescales and compute
                % change since last computation
                Ds = [Ds {currentParams.D}];
                gams = [gams {currentParams.gamma}];
            end
    end
    
    
    % ======================
    % Mean parameter, d
    % ======================
    
    currentParams.d.cov = 1./(prior.d.beta + NT*currentParams.phi.mean);
    currentParams.d.mean = diag(currentParams.d.cov) * diag(currentParams.phi.mean) ...
        * sum((Y - blkdiag(currentParams.C.means{:}) * X),2);
    dd = diag(diag(currentParams.d.cov) ...
        + currentParams.d.mean*currentParams.d.mean'); % Second moment of d
    % Update zero-centered observations
    Y0 = Y - repmat(currentParams.d.mean,1,NT);
    
    % ======================
    % Loading matrices, C
    % ======================
    logdetC = 0; % To be used in calculation of lower bound
    for groupIdx = 1:numGroups
        currGroup = block_idxs{groupIdx};
        phi_m = currentParams.phi.mean(currGroup(1):currGroup(2));
        alpha_m = diag(currentParams.alpha.mean(groupIdx,:));  % (xDim x xDim) array
        XY_m = Xs{groupIdx} * Y0(currGroup(1):currGroup(2),:)';
        for yIdx = 1:yDims(groupIdx)
            % Covariance
            covC = inv(alpha_m + phi_m(yIdx) .* XX{groupIdx});
            covC = 0.5 * (covC + covC'); % Ensure symmetry
            logdetC = logdetC + logdet(covC);
            if saveCcov
                currentParams.C.covs{groupIdx}{yIdx} = covC;
            end
            % Mean
            currentParams.C.means{groupIdx}(yIdx,:) ...
                = (phi_m(yIdx) * covC * XY_m(:,yIdx))';
            % Second moment
            currentParams.C.moments{groupIdx}{yIdx} = covC ...
                + currentParams.C.means{groupIdx}(yIdx,:)' * currentParams.C.means{groupIdx}(yIdx,:);
        end
    end
    
    % ======================
    % ARD parameters, alpha
    % ======================
    
    for groupIdx = 1:numGroups
        CC_m = zeros(xDim);  % Second moments of C_m (m = groupIdx)
        for yIdx = 1:yDims(groupIdx)
            CC_m = CC_m + currentParams.C.moments{groupIdx}{yIdx}; 
        end
        currentParams.alpha.b(groupIdx,:) = (prior.alpha.b + diag(CC_m)./2)';
        currentParams.alpha.mean(groupIdx,:) ...
            = currentParams.alpha.a(groupIdx) ./ currentParams.alpha.b(groupIdx,:); 
    end
    if (rem(i, freqParam) == 0) || (i == maxIters)
        alphas = [alphas {currentParams.alpha.mean}];
    end
    
    % ======================
    % Noise precisions, phi
    % ======================
    
    for groupIdx = 1:numGroups
        currGroup = block_idxs{groupIdx};
        phi_bm = zeros(yDims(groupIdx),1);
        dd_m = dd(currGroup(1):currGroup(2));
        d_m = currentParams.d.mean(currGroup(1):currGroup(2));
        Y_m = Y(currGroup(1):currGroup(2),:);
        Y0_m = Y0(currGroup(1):currGroup(2),:);
        Y2_m = Y2(currGroup(1):currGroup(2),:);
        for yIdx = 1:yDims(groupIdx)
            phi_bm(yIdx) = prior.phi.b + 0.5 * ( ...
                           NT * dd_m(yIdx) + sum(Y2_m(yIdx,:) - 2 * d_m(yIdx) * Y_m(yIdx,:),2) ...
                           - 2.*currentParams.C.means{groupIdx}(yIdx,:) * Xs{groupIdx} * Y0_m(yIdx,:)' ...
                           + trace(currentParams.C.moments{groupIdx}{yIdx} * XX{groupIdx}) ...
                           );
        end
        currentParams.phi.b(currGroup(1):currGroup(2)) = phi_bm;
        currentParams.phi.mean(currGroup(1):currGroup(2)) ...
            = currentParams.phi.a ./ phi_bm;
    end
    % Set minimum private variance
    currentParams.phi.mean = min(1./varFloor, currentParams.phi.mean);
    currentParams.phi.b = currentParams.phi.a ./ currentParams.phi.mean;
    
    % ======================
    % Compute lower bound
    % ======================
    
    if ~isnan(lbi)
        lbold = lbi;
    end
    
    if getLB
        % Likelihood term
        log_phi = digammaa_phi - log(currentParams.phi.b); % (yDim x 1) array
        lbi = const_lik + (NT/2) * sum(log_phi) ...
            - sum(currentParams.phi.mean .* (currentParams.phi.b - prior.phi.b));

        % X KL term
        lbi = lbi + (numGroups*xDim*NT)/2 + res.lb_gp + 0.5*lbterms.logdet_SigX;

        log_alpha = nan(numGroups,xDim);
        for groupIdx = 1:numGroups
            log_alpha(groupIdx,:) ...
                = digammaa_alpha(groupIdx) - log(currentParams.alpha.b(groupIdx,:));        
        end

        % C KL term
        lbi = lbi + 0.5*logdetC;
        for groupIdx = 1:numGroups
            lbi = lbi + (yDims(groupIdx)/2)*sum(log_alpha(groupIdx,:));
            for yIdx = 1:yDims(groupIdx)
                lbi = lbi + 0.5*trace(eye(xDim) ...
                    - currentParams.C.moments{groupIdx}{yIdx}*diag(currentParams.alpha.mean(groupIdx,:))); 
            end
        end

        % alpha KL term
        lbi = lbi + (numGroups*xDim) * (alogb_alpha - loggammaa_alpha_prior);
        for groupIdx = 1:numGroups
            lbi = lbi + sum( ...
                -currentParams.alpha.a(groupIdx) .* log(currentParams.alpha.b(groupIdx,:)) ...
                - prior.alpha.b .* currentParams.alpha.mean(groupIdx,:) ...
                + (prior.alpha.a - currentParams.alpha.a(groupIdx)) .* (log_alpha(groupIdx,:)) ...
                ) ...
                + xDim * loggammaa_alpha_post(groupIdx) ...
                + xDim * currentParams.alpha.a(groupIdx);
        end

        % phi KL term
        lbi = lbi + yDim * (alogb_phi + loggammaa_phi_post - loggammaa_phi_prior + currentParams.phi.a) ...
            + sum( ...
                  -currentParams.phi.a .* log(currentParams.phi.b) ...
                  - prior.phi.b .* currentParams.phi.mean ...
                  + (prior.phi.a - currentParams.phi.a) .* (digammaa_phi - log(currentParams.phi.b)) ...
              );

        % d KL term
        lbi = lbi + yDim/2 + (yDim/2)*log(prior.d.beta) + 0.5*logdet(diag(currentParams.d.cov)) ...
            - 0.5 * prior.d.beta * sum(dd);
    else
        lbi = nan;    
    end
    
    % Finish tracking EM iteration time
    tEnd    = toc;
    iterTime = [iterTime tEnd];
    
    % Check stopping conditions or errors
    if verbose
        if getLB
            fprintf('EM iteration %3d of %d        lb %f\r', i, maxIters, lbi);
        else
            fprintf('EM iteration %3d of %d\r', i, maxIters);
        end
    end
    
    lb = [lb lbi];
    
    if i <= 2
        lbbase = lbi;
    elseif (lbi < lbold)
        flags.DecreasingLowerBound = 1;
        disp('Error: Decreasing lower bound');
    elseif ((lbi - lbbase) < (1+tol)*(lbold - lbbase))
        break;
    end
    
end

if (length(lb) < maxIters) && (xDim > 0)
    flags.Convergence = 1;
end

if verbose
    if flags.Convergence == 1
        fprintf('Lower bound converged after %d iterations.\n', length(lb));
    elseif (length(lb) < maxIters) && (xDim <= 0)
        fprintf('Fitting stopped because no significant latent dimensions remain.\n');
    else
        fprintf('Fitting stopped after maxIters (%d) was reached.\n', maxIters);
    end 
end

if any(currentParams.phi.mean == 1./varFloor)
    flags.PrivateVarianceFloor = 1; 
end

% Save tracked parameters
trackedParams.lb = lb;
trackedParams.iterTime = iterTime;
trackedParams.Ds = Ds;
trackedParams.gams = gams;
trackedParams.alphas = alphas;

% Remove estimates of X.cov, if desired
if ~saveXcov
    currentParams = rmfield(currentParams,'X');
end
