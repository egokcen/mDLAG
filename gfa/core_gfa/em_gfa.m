function outparams = em_gfa(Ys, xDim, varargin)
%
% outparams = em_gfa(Ys, xDim, ...)
%
% Description: Fit a group factor analysis (GFA) model using a variational
%              EM algorithm with mean-field approximation.
%
% Arguments:
%
%     Required:
%
%     Ys          -- (1 x numGroups) cell array; list of data matrices 
%                    {(y1Dim x N), (y2Dim x N), ...}
%     xDim        -- int; number of factors
%
%     Optional:
%
%     R        -- string or int; 'full' for full-rank alpha model, a
%                 positive integer giving the rank of alpha otherwise.
%     include_mu -- logical; set true to include mean parameters in
%                   low-rank alpha model. (default: true)
%     prior    -- structure with the following fields:
%                   d.beta  -- positive float; precision of mean parameter
%                              generative model (Gaussian)
%                   phi.a   -- positive float; 'a' shape parameter of
%                              observation precision (phi) generative 
%                              model (Gamma with mean a/b)
%                   phi.b   -- positive float; 'b' scale parameter of
%                              observation precision (phi) generative
%                              model (Gamma with mean a/b)
%                   alpha.a -- positive float; 'a' shape parameter of 
%                              alpha parameter generative model 
%                              (Gamma with mean a/b)
%                   alpha.b -- positive float; 'b' scale parameter of
%                              alpha parameter generative model 
%                              (Gamma with mean a/b)
%                   uv.lbda -- positive float; regularization parameter
%                              for low-rank alpha model
%                   (default: uv.lbda, 0.1; all other values, 1e-12)
%
%     tol      -- float; stopping criterion for EM (default: 1e-8)
%     maxIters -- int; maximum number of EM iterations (default: 1e8)
%     verbose  -- boolean; specifies whether to display status messages
%                 (default: false)
%     randomSeed -- int or string; seed the random number generator, 
%                   for reproducible initializations (default: 'shuffle')
%     minVarFrac -- float; fraction of overall data variance for each 
%                   observed dimension to set as the private variance 
%                   floor. (default: 0.001)
%     pruneX   -- logical; set true to remove dimensions from X that
%                 become inactive. Can speed up EM runtime and improve
%                 memory efficiency. (default: true)
%     saveX    -- logical; set true to save posterior estimates of latent
%                 variables X. For large datasets, X.mean may be very
%                 large. (default: false)
%     saveCcov -- logical; set true to save posterior covariance and 
%                 of C. For large yDim and xDim, these structures can use
%                 a lot of memory. (default: false)
%     saveFitProgress -- logical; set true to save lower bound and
%                        iteration time each EM iteration. (default: true)
%
% Outputs:
%     outparams -- structure containing values relevant to the fitted
%                  model:
%         xDim       -- int; number of factors
%         yDims      -- (1 x numGroups) array; dimensionalities of each 
%                       observed group
%         X.mean     -- (xDim x N) array; posterior mean of latent
%                       variables
%         X.cov      -- (xDim x xDim) array; posterior covariance of latent
%                       variables
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
%         if R == 'full' (i.e., full-rank alpha model)
%             alpha.a    -- (numGroups x 1) array; shape parameters of 
%                           alpha posterior
%             alpha.b    -- (numGroups x xDim) array; scale parameters of 
%                           alpha posterior
%             alpha.mean -- (numGroups x xDim) array; mean precisions of
%                           loading weights (for ARD); alpha.a ./ alpha.b
%         else           (i.e., low-rank alpha model)
%             alpha.mean -- (numGroups x xDim) array; precisions of
%                           loading weights (for ARD)
%           if R > 0
%             alpha.U    -- (numGroups x R) array; low-rank factorization
%                           of alpha
%             alpha.V    -- (xDim x R) array; low-rank factorization of
%                           alpha
%           if include_mu
%             alpha.mu_u -- (numGroups x 1) array; mean offset in low-rank
%                           factorization
%             alpha.mu_v -- (xDim x 1) array; mean offset in low-rank
%                           factorization
%         phi.a      -- float; shape parameter of phi posterior
%         phi.b      -- (yDim x 1) array; scale parameters of phi posterior
%         phi.mean   -- (yDim x 1) array; mean precisions of observations; 
%                       alpha.a ./ alpha.b
%         lb         -- (1 x numIters) array; variational lower bound at
%                       each iteration
%         iterTime   -- (1 x numIters) array; computation time for each EM
%                       iteration
%         flags      -- structure containing various status messages about
%                       the fitting procedure:
%             Convergence          -- logical; 1 if lower bound converged
%                                     before reaching maxIters EM
%                                     iterations; 0 otherwise
%             DecreasingLowerBound -- logical; 1 if lower bound decreased 
%                                     during fitting; 0 otherwise.
%             PrivateVarianceFloor -- logical; 1 if private variance
%                                     floor was used on any values of phi;
%                                     0 otherwise.
%             xDimsRemoved         -- int; Number of latent dimensions
%                                     removed (if pruneX is true) due to 
%                                     low variance.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     06 Jun 2022 -- Initial full revision.
%     01 Jul 2022 -- Added low-rank alpha model. Cleaned up comments,
%                    documentation.
%     09 Jul 2022 -- Added option to save posterior estimates of X.
%                    Added R = 0 compatibility.
%     12 Jul 2022 -- Added option to remove mean parameters from low-rank
%                    alpha model.
%     13 Jul 2022 -- Added variance floor option.
%     23 Oct 2022 -- Added the performance-enhancing pruneX and saveCcov 
%                    options.

R = 'full';
include_mu = true;
prior.d.beta = 1e-12;
prior.phi.a = 1e-12;
prior.phi.b = 1e-12;
prior.alpha.a = 1e-12;
prior.alpha.b = 1e-12;
prior.uv.lbda = 0.1;
tol = 1e-8; 
maxIters = 1e8;
verbose = false;
randomSeed = 'shuffle';
minVarFrac = 0.001;
pruneX = true;
saveX = false;
saveCcov = false;
saveFitProgress = true;
assignopts(who, varargin);

% Large alpha leads to this warning, but doesn't seem to cause accuracy
% issues in practice.
warning('off', 'MATLAB:nearlySingularMatrix');

rng(randomSeed);

fname = 'grad_uv'; % Name of function to evaluate U, V gradients, if relevant
MAXGRADITERS = -100; % Maximum number of gradient iterations in minimize.m

numGroups = length(Ys);
% Get individual group sizes
yDims = zeros(1,numGroups);
for groupIdx = 1:numGroups
    yDims(groupIdx) = size(Ys{groupIdx},1); 
end

% Useful for extracting correct-sized blocks from matrices later
block_idxs = get_block_idxs(yDims);

if ~isequal(R, 'full')
    
    if R <= 0 && ~include_mu
        fprintf('Error: If R = 0, then mu must be included in the low-rank alpha model.\n');
        outparams = [];
        return
    end

    if R >= min([xDim numGroups])
        if verbose
            fprintf('Specified R corresponds to full rank solution.\n');
        end
        R = 'full';
    end
end

% Construct joint data matrix
Y = cat(1,Ys{:});
[yDim, N] = size(Y);
Y2 = Y.^2;

% Get total variance of each group
covY = cov(Y', 1);
totalvars = nan(1,numGroups);
for groupIdx = 1:numGroups
    currGroup = block_idxs{groupIdx};
    totalvars(groupIdx) = trace(covY(currGroup(1):currGroup(2),currGroup(1):currGroup(2)));
end
varFloor = minVarFrac * diag(covY);

% Define constants
% Latent magnitude must be greater than this value to remain in the model
PRUNEX_TOL = 1e-7;

% Constant factors in elbo
const = -((yDim*N)/2) * log(2*pi);

% Initialize the model

% Latent variables
X.mean = randn(xDim,N);         % Mean
X.cov = eye(xDim);              % Covariance

% Mean parameter
d.mean = mean(Y,2);                      % Mean
d.cov = (1/prior.d.beta).*ones(yDim,1);  % Covariance

% Noise precisions
phi.a = prior.phi.a + N/2;
phi.b = prior.phi.b .* ones(yDim,1);
alogb_phi = prior.phi.a * log(prior.phi.b);
loggammaa_phi_prior = gammaln(prior.phi.a);
loggammaa_phi_post = gammaln(phi.a);
digammaa_phi = psi(phi.a);
phi.mean = 1./diag(covY);

% Loading matrices
% No need for random initialization, since C will be the first update below
C.means = cell(1,numGroups);
if saveCcov
    C.covs = cell(1,numGroups);
end
CC = cell(1,numGroups);
for groupIdx = 1:numGroups
    C.means{groupIdx} = zeros(yDims(groupIdx),xDim);
    covY_m = cov(Ys{groupIdx}', 1);
    if rank(covY_m) == yDims(groupIdx)
        scale = exp(2*sum(log(diag(chol(covY_m))))/yDim); % Scale by determinant of covY_m
    else
        % covY_m may not be full rank because N < yDim
        if verbose
            fprintf('WARNING: Data matrix for group %d is not full rank.\n', groupIdx);
        end
        r     = rank(covY_m);
        e     = sort(eig(covY_m), 'descend');
        scale = geomean(e(1:r));
    end
    C.means{groupIdx} = randn(yDims(groupIdx),xDim)*sqrt(scale/xDim);
    if saveCcov
        C.covs{groupIdx} = cell(yDims(groupIdx),1);
    end
    CC{groupIdx} = cell(yDims(groupIdx),1);
    for yIdx = 1:yDims(groupIdx)
        covC = zeros(xDim);
        if saveCcov
            C.covs{groupIdx}{yIdx} = covC;
        end
        CC{groupIdx}{yIdx} = covC ...
            + C.means{groupIdx}(yIdx,:)'*C.means{groupIdx}(yIdx,:);
    end
end
Call = vertcat(C.means{:});

% ARD parameters
if isequal(R, 'full')
    alpha.a = prior.alpha.a + (yDims./2)';  % (numGroups x 1) array
    alpha.b = prior.alpha.b .* ones(numGroups,xDim);
    alogb_alpha = prior.alpha.a * log(prior.alpha.b);
    loggammaa_alpha_prior = gammaln(prior.alpha.a);
    loggammaa_alpha_post = gammaln(alpha.a); % (numGroups x 1) array
    digammaa_alpha = psi(alpha.a);          % (numGroups x 1) array
end
alpha.mean = nan(numGroups,xDim);
% Scale ARD parameters to match the data
for groupIdx = 1:numGroups
    CC_m = zeros(xDim);
    for yIdx = 1:yDims(groupIdx)
        CC_m = CC_m + CC{groupIdx}{yIdx}; 
    end
    alpha.mean(groupIdx,:) = (yDims(groupIdx)./diag(CC_m))';
end

% Low-rank alpha factorization
if ~isequal(R, 'full')
    uv = [];
    if R > 0
        alpha.U = abs(randn(numGroups,R));
        alpha.V = abs(randn(xDim,R));
        % Vectorized form of UV parameters for optimization
        uv = [uv; alpha.U(:); alpha.V(:)];
    end
    
    if include_mu
        alpha.mu_u = zeros(numGroups,1);
        alpha.mu_v = zeros(xDim,1);
        % Vectorized form of UV parameters for optimization
        uv = [uv; alpha.mu_u; alpha.mu_v];
    end
    lenU = numGroups*R;
    lenV = xDim*R;
    C2 = nan(numGroups,xDim); % Second moments of C
end

% Variational lower bound
lbi = -Inf;      % Initial lower bound
lb = [];         % Lower bound at each iteration
iterTime = [];   % Time it takes to complete each iteration

% Initialize status flags
outparams.flags.Convergence = 0;
outparams.flags.DecreasingLowerBound = 0;
outparams.flags.PrivateVarianceFloor = 0;
outparams.flags.xDimsRemoved = 0;

% Begin EM iterations
for i = 1:maxIters
    
    tic; % For tracking the computation time of each iteration
    
    % Check if any latent states need to be removed
    if pruneX
        keptXdims = find(mean(X.mean.^2,2) > PRUNEX_TOL);
        if length(keptXdims) < xDim
            % Remove insignificant latent dimensions from all relevant
            % model parameters.
            
            % ARD parameters
            alpha.mean = alpha.mean(:,keptXdims);
            if isequal(R, 'full')
                alpha.b = alpha.b(:,keptXdims);
            else
                uv = [];
                if R > 0
                    alpha.V = alpha.V(keptXdims,:);
                    % Vectorized form of UV parameters for optimization
                    uv = [uv; alpha.U(:); alpha.V(:)];
                end

                if include_mu
                    alpha.mu_v = alpha.mu_v(keptXdims);
                    % Vectorized form of UV parameters for optimization
                    uv = [uv; alpha.mu_u; alpha.mu_v];
                end
                lenV = length(keptXdims)*R;
                C2 = C2(:,keptXdims); % Second moments of C
            end
            
            % Loading matrix, C
            for groupIdx = 1:numGroups
                C.means{groupIdx} = C.means{groupIdx}(:,keptXdims);
                for yIdx = 1:yDims(groupIdx)
                    if isfield(C, 'covs')
                        C.covs{groupIdx}{yIdx} ...
                            = C.covs{groupIdx}{yIdx}(keptXdims,keptXdims);
                    end
                    CC{groupIdx}{yIdx} ...
                        = CC{groupIdx}{yIdx}(keptXdims,keptXdims);
                end
            end
            Call = vertcat(C.means{:});
            
            % Latent variables, X
            X.mean = X.mean(keptXdims,:);
            X.cov = X.cov(keptXdims,keptXdims);
            
            outparams.flags.xDimsRemoved = outparams.flags.xDimsRemoved ...
                + (xDim - length(keptXdims));
            xDim = length(keptXdims);
            if xDim < 0
                % Stop fitting if no significant latent dimensions remain.
                break;
            end
        end
    end
    
    % Mean parameter, d
    d.cov = 1./(prior.d.beta + N*phi.mean);
    d.mean = diag(d.cov) * diag(phi.mean) * sum((Y - Call * X.mean),2);
    dd = diag(diag(d.cov) + d.mean*d.mean');
    
    % Latent variables, X
    
    % Zero-center observations with current estimate of means
    Y0 = Y - repmat(d.mean,1,N);
    
    % Covariance
    X.cov = eye(xDim);
    for groupIdx = 1:numGroups
        currGroup = block_idxs{groupIdx};
        phi_m = phi.mean(currGroup(1):currGroup(2));
        for yIdx = 1:yDims(groupIdx)
            X.cov = X.cov + phi_m(yIdx).*CC{groupIdx}{yIdx};
        end
    end
    X.cov = inv(X.cov);
    X.cov = 0.5 * (X.cov + X.cov'); % Ensure symmetry
    % Mean
    X.mean = X.cov * Call' * diag(phi.mean) * Y0;
    % Second moment
    XX = N.*X.cov + X.mean*X.mean';
   
    % Loading matrices, C
    
    % Intermediate computations
    
    % Correlation matrix between current estimate of latent variables and
    % zero-centered observations.
    XY = X.mean*Y0';   % (xDim x yDim) array
    logdetC = 0;       % To be used in calculation of lower bound
    for groupIdx = 1:numGroups
        currGroup = block_idxs{groupIdx};
        phi_m = phi.mean(currGroup(1):currGroup(2));
        alpha_m = diag(alpha.mean(groupIdx,:));  % (xDim x xDim) array
        XY_m = XY(:,currGroup(1):currGroup(2));
        for yIdx = 1:yDims(groupIdx)
            % Covariance
            covC = inv(alpha_m + phi_m(yIdx) .* XX);
            covC = 0.5 * (covC + covC'); % Ensure symmetry
            logdetC = logdetC + logdet(covC);
            if saveCcov
                C.covs{groupIdx}{yIdx} = covC;
            end
            % Mean
            C.means{groupIdx}(yIdx,:) = (phi_m(yIdx) * covC * XY_m(:,yIdx))';
            % Second moment
            CC{groupIdx}{yIdx} = covC ...
                + C.means{groupIdx}(yIdx,:)' * C.means{groupIdx}(yIdx,:);
        end
    end
    Call = vertcat(C.means{:});
    
    % ARD parameters, alpha
    if isequal(R, 'full')
        for groupIdx = 1:numGroups
            CC_m = zeros(xDim);
            for yIdx = 1:yDims(groupIdx)
                CC_m = CC_m + CC{groupIdx}{yIdx}; 
            end
            alpha.b(groupIdx,:) = (prior.alpha.b + diag(CC_m)./2)';
            alpha.mean(groupIdx,:) = alpha.a(groupIdx) ./ alpha.b(groupIdx,:); 
        end
    else
        for groupIdx = 1:numGroups
            CC_m = zeros(xDim);
            for yIdx = 1:yDims(groupIdx)
                CC_m = CC_m + CC{groupIdx}{yIdx}; 
            end
            C2(groupIdx,:) = diag(CC_m);
        end
        
        % Minimize the lower bound with respect to U, V, mu_u, mu_v
        [uv, lb_uv, ~] = minimize(uv, fname, MAXGRADITERS, C2, yDims, R, prior.uv.lbda, include_mu);
        
        % Unpack output of minimization
        lb_uv = -lb_uv(end);
        alpha.mean = zeros(numGroups,xDim);
        if R > 0
            alpha.U = reshape(uv(1:lenU), size(alpha.U));
            alpha.V = reshape(uv(lenU+1:lenU+lenV), size(alpha.V));
            alpha.mean = alpha.mean + alpha.U * alpha.V';
        end
        if include_mu
            alpha.mu_u = uv(lenU+lenV+1:lenU+lenV+numGroups);
            alpha.mu_v = uv(lenU+lenV+numGroups+1:end);
            alpha.mean = alpha.mean + alpha.mu_u*ones(1,xDim) + ones(numGroups,1)*alpha.mu_v';
        end
        alpha.mean = exp(alpha.mean);
    end
    
    % Noise precisions, phi
    for groupIdx = 1:numGroups
        currGroup = block_idxs{groupIdx};
        phi_bm = zeros(yDims(groupIdx),1);
        dd_m = dd(currGroup(1):currGroup(2));
        d_m = d.mean(currGroup(1):currGroup(2));
        Y_m = Y(currGroup(1):currGroup(2),:);
        Y0_m = Y0(currGroup(1):currGroup(2),:);
        Y2_m = Y2(currGroup(1):currGroup(2),:);
        for yIdx = 1:yDims(groupIdx)
            phi_bm(yIdx) = prior.phi.b + 0.5 * ( ...
                           N * dd_m(yIdx) + sum(Y2_m(yIdx,:) - 2 * d_m(yIdx) * Y_m(yIdx,:),2) ...
                           - 2.*C.means{groupIdx}(yIdx,:) * X.mean * Y0_m(yIdx,:)' ...
                           + trace(CC{groupIdx}{yIdx} * XX) ...
                           );
        end
        phi.b(currGroup(1):currGroup(2)) = phi_bm;
        phi.mean(currGroup(1):currGroup(2)) = phi.a ./ phi_bm;
    end
    % Set minimum private variance
    phi.mean = min(1./varFloor, phi.mean);
    phi.b = phi.a ./ phi.mean;
    
    % Compute lower bound
    lbold = lbi;
    
    % Likelihood term
    log_phi = digammaa_phi - log(phi.b); % (yDim x 1) array
    lbi = const + (N/2) * sum(log_phi) - sum(phi.mean .* (phi.b - prior.phi.b));
    
    % X KL term
    lbi = lbi + (N/2)*logdet(X.cov) + (1/2)*trace(N.*eye(xDim) - XX);
    
    if isequal(R, 'full')
        
        log_alpha = nan(numGroups,xDim);
        for groupIdx = 1:numGroups
            log_alpha(groupIdx,:) = digammaa_alpha(groupIdx) - log(alpha.b(groupIdx,:));        
        end
        
        % C KL term
        lbi = lbi + 0.5*logdetC;
        for groupIdx = 1:numGroups
            lbi = lbi + (yDims(groupIdx)/2)*sum(log_alpha(groupIdx,:));
            for yIdx = 1:yDims(groupIdx)
                lbi = lbi + 0.5*trace(eye(xDim) ...
                    - CC{groupIdx}{yIdx}*diag(alpha.mean(groupIdx,:))); 
            end
        end
    
        % alpha KL term
        lbi = lbi + (numGroups*xDim) * (alogb_alpha - loggammaa_alpha_prior);
        for groupIdx = 1:numGroups
            lbi = lbi + sum( ...
                -alpha.a(groupIdx) .* log(alpha.b(groupIdx,:)) ...
                - prior.alpha.b .* alpha.mean(groupIdx,:) ...
                + (prior.alpha.a - alpha.a(groupIdx)) .* (log_alpha(groupIdx,:)) ...
                ) ...
                + xDim * loggammaa_alpha_post(groupIdx) ...
                + xDim * alpha.a(groupIdx);
        end
        
    else
        
        % C KL term
        lbi = lbi + 0.5*yDim*xDim + 0.5*logdetC;
        lbi = lbi + lb_uv;
    end
    
    % phi KL term
    lbi = lbi + yDim * (alogb_phi + loggammaa_phi_post - loggammaa_phi_prior + phi.a) ...
        + sum( ...
              -phi.a .* log(phi.b) ...
              - prior.phi.b .* phi.mean ...
              + (prior.phi.a - phi.a) .* (digammaa_phi - log(phi.b)) ...
          );
    
    % d KL term
    lbi = lbi + yDim/2 + (yDim/2)*log(prior.d.beta) + 0.5*logdet(diag(d.cov)) ...
        - 0.5 * prior.d.beta * sum(dd);
    
    % Finish tracking EM iteration time
    tEnd    = toc;
    iterTime = [iterTime tEnd];
    
    % Check stopping conditions or errors
    if verbose
        fprintf('EM iteration %3d of %d        lb %f\r', i, maxIters, lbi);
    end
    
    lb = [lb lbi];
    
    if i <= 2
        lbbase = lbi;
    elseif (lbi < lbold)
        outparams.flags.DecreasingLowerBound = 1;
        disp('Error: Decreasing lower bound');
    elseif ((lbi - lbbase) < (1+tol)*(lbold - lbbase))
        break;
    end
    
end

if (length(lb) < maxIters) && (xDim > 0)
    outparams.flags.Convergence = 1;
end

if verbose
    if outparams.flags.Convergence == 1
        fprintf('Lower bound converged after %d iterations.\n', length(lb));
    elseif (length(lb) < maxIters) && (xDim <= 0)
        fprintf('Fitting stopped because no significant latent dimensions remain.\n');
    else
        fprintf('Fitting stopped after maxIters (%d) was reached.\n', maxIters);
    end 
end

% Collect outputs
outparams.xDim = xDim;
outparams.yDims = yDims;
if saveX
    outparams.X.mean = X.mean;
    outparams.X.cov = X.cov;
end
outparams.d.mean = d.mean;
outparams.d.cov = d.cov;
outparams.C.means = C.means;
if saveCcov
    outparams.C.covs = C.covs;
end
outparams.C.moments = CC;
outparams.alpha.mean = alpha.mean;
if isequal(R, 'full')
    outparams.alpha.a = alpha.a;
    outparams.alpha.b = alpha.b;
else
    if R > 0
        outparams.alpha.U = alpha.U;
        outparams.alpha.V = alpha.V;
    end
    if include_mu
        outparams.alpha.mu_u = alpha.mu_u;
        outparams.alpha.mu_v = alpha.mu_v;
    end
end
outparams.phi.mean = phi.mean;
outparams.phi.a = phi.a;
outparams.phi.b = phi.b;
if saveFitProgress
    outparams.lb = lb;
    outparams.iterTime = iterTime;
end
if any(outparams.phi.mean == 1./varFloor)
    outparams.flags.PrivateVarianceFloor = 1; 
end
