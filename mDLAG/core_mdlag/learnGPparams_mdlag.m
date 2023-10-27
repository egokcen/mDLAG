function res = learnGPparams_mdlag(seq, params, constraints, varargin)
%
% res = learnGPparams_mdlag(seq, params, constraints, ...)
%
% Description: Update parameters of GP state model given inferred
%              latent states.
%              NOTE: Learning GP noise variance (eps) is currently
%                    unsupported.
%
% Arguments:
%
%     Required:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
%                     xsm          -- ((numGroups*xDim) x T) array; 
%                                     posterior mean at each timepoint
% 
%     params -- Structure containing mDLAG model parameters.
%               Contains the fields
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
%     constraints -- Structure containing constraints on GP parameters:
%                      maxDelay -- float; maximum delay magnitude, in units
%                                  of time steps
%                      minGamma -- float; minimum gamma value 
%                                  (binWidth^2/tau^2), unitless
%
%     Optional:
%
%     MAXITERS -- int; maximum number of line searches (if >0), 
%                 maximum number of function evaluations (if <0), 
%                 for minimize.m (default:-10)
%     verbose  -- logical that specifies whether to display status messages
%                 (default: false)
%
% Outputs:
%
%     res -- Structure containing the updated GP state model parameters:
%            gamma, D. Also includes the number of iterations and
%            value of the cost function after updating these parameters via
%            gradient descent.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     19 Oct 2022 -- Initial full revision. 

MAXITERS  = -10; % for minimize.m
verbose   = false;
assignopts(who, varargin);

xDim         = params.xDim;
yDims        = params.yDims;
numGroups    = length(yDims);

switch params.covType
    case 'rbf'
        % If there's more than one type of parameter, put them in the
        % second row of oldParams.
        oldParams = [params.gamma; params.D];
        fname     = 'grad_GPparams';
    % Optional: Insert other covariance functions here    
end

precomp = makePrecomp_GP(seq,params);

% Loop once for each state dimension (each GP)
gamma = zeros(1,xDim);
D = zeros(numGroups, xDim);
res.lb_gp = 0; % Value of lower bound component involving GP params
for i = 1:xDim
    const = [];
    
    % We're not learning GP noise variance, for now
    const.eps = params.eps(i);
    const.minGamma = constraints.minGamma;
    const.maxDelay = constraints.maxDelay;
    
    switch fname                
        case 'grad_GPparams'
            % Change of variables for constrained optimization
            init_gam = log(oldParams(1,i) - constraints.minGamma);
            % We don't include delays to the first group in the optimization
            init_delay = reshape(oldParams(3:end,i), numGroups-1, 1);
            init_delay = log(constraints.maxDelay + init_delay) ...
                       - log(constraints.maxDelay - init_delay);
            init_p = [init_gam; init_delay];
    end   
    
    % This does the heavy lifting
    [res_p, f_gp, grad_iters] =...
        minimize(init_p, fname, MAXITERS, precomp(i), const);
    
    switch params.covType
        case 'rbf'
            switch fname                
                case 'grad_GPparams'
                    gamma(i) = exp(res_p(1)) + constraints.minGamma;
                    D(2:end,i) = constraints.maxDelay.*tanh(res_p(2:end)./2);
                    res.lb_gp = res.lb_gp + -f_gp(end);
            end        
    end    
    
    if verbose
        fprintf('\nConverged p; xDim:%d, p:%s', i, mat2str(res_p, 3));
    end
end

res.D = D;
res.gamma = gamma;

end

