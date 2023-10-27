function seq = generate_latents_mdlag(params, T, N, varargin)
% 
% seq = generate_latents_mdlag(params, T, N,...)
%
% Description: Generate N independent sequences of length T samples, 
%              according to a zero-mean Gaussian Process defined by the 
%              mDLAG state model.
%
% Arguments:
%
%     Required:
%
%     params  -- Structure containing mDLAG model parameters.
%                Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    Cs      -- (1 x numGroups) cell array; List of factor 
%                               loadings 
%                               {(y1Dim x xDim), (y2Dim x xDim), ...}
%                    alphas  -- (numGroups x xDim) array; alpha parameter 
%                               values
%                    phis    -- (1 x numGroups) cell array; List of 
%                               observation precision parameters 
%                               {(y1Dim x 1), (y2Dim x 1), ...}
%                    ds      -- (1 x numGroups) cell array; List of data 
%                               means
%                               {(y1Dim x 1), (y2Dim x 1), ...}
%                    gamma   -- (1 x xDim) array; GP timescales in units of 
%                               time are given by 'binWidth ./ sqrt(gamma)'                                                    
%                    eps     -- (1 x xDim) array; GP noise variances
%                    D       -- (numGroups x xDim) array; delays from 
%                               latents to observed variables. NOTE: Delays
%                               are reported as (real-valued) number of 
%                               time-steps.
%                    xDim    -- int; number of latent variables
%                    yDims   -- (1 x numGroups) array; dimensionalities of 
%                               each observed group
%     T         -- (1 x N) int array; T(n) gives the number of samples 
%                  (time length) for sequence n. If T is a scalar
%                  (length-1), then all sequences will be length-T.
%     N         -- int; number of sequences
%
%     Optional:
%
%     latentfield -- string; Name of data field in seq (default: 'xsm')
%     verbose     -- logical; Print status info (default: false)
%
% Outputs:
%     seq -- structure whose nth entry (corresponding to the nth 
%            sequence) has fields
%                trialId   -- int; unique trial (sequence) identifier  
%                T -- int; number of timesteps
%                (latentfield) -- (numGroups*xDim x T) array; delayed
%                                 latent sequences
%
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     27 Sep 2022 -- Initial full revision.

latentfield = 'xsm';
verbose = false;
extraOpts = assignopts(who, varargin);
   
% If T is a scalar, then give all sequences the same length
if length(T) <= 1
    T = repmat(T,1,N);
end

% Group trials of same length together
Tu = unique(T);

% Initialize output structure
numGroups = length(params.yDims);
xDim = params.xDim;
for n = 1:N
    seq(n).trialId = n;
    seq(n).T = T(n);
    seq(n).(latentfield) = nan(numGroups*xDim,T(n));
end
    
% Generate all trials of the same length
for j = 1:length(Tu)
    Tj = Tu(j);
    if verbose
        fprintf('Generating all trials of length T = %d..\n', Tj);
    end
    nList = find(T == Tj);
    K_big = construct_K_mdlag(params, Tj); % GP kernel matrix
    for n = nList
        if verbose
            fprintf('    Trial n = %d...\n', nList(n));
        end
        seq(n).(latentfield) = reshape(mvnrnd(zeros(1,size(K_big,1)), K_big),[],Tj);
    end
end
