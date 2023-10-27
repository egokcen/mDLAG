function precomp = makePrecomp_GP(seq,params)
%
% precomp = makePrecomp_GP(seq,params)
%
% Description: Precompute posterior second moment for the GP parameter
%              update.
%
% Arguments:
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
%
% Outputs
%     precomp -- Structure whose ith entry contains precomputations for 
%                the i-th latent state:
%                Tdif -- (numGroups*T x numGroups*T) array; matrix of time
%                        differences.
%                Tall -- (1 x N) array; List of trial lengths
%                params -- structure of the same format as 'params' above.
%                Tu   -- structure whose jth entry, corresponding to a
%                        group of trials of the same length, contains 
%                        the following:
%                        nList -- List of trial IDs belonging to this group
%                        T     -- int; Length of all trials in this group
%                        numTrials -- int; Number of trials in this group
%                        Xmoment -- (numGroups*T x numGroups*T) array;
%                                   Precomputed posterior second moment of X 
%                                   based on the observations in this group            
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     20 Oct 2022 -- Initial full revision. 

% Initialize relevant variables
xDim = params.xDim;
yDims = params.yDims;
numGroups = length(yDims);
Tall = [seq.T];
Tmax = max(Tall);
Tdif = repmat(1:Tmax,numGroups,1);
Tdif = repmat(Tdif(:)',numGroups*Tmax,1) - repmat(Tdif(:),1,numGroups*Tmax);

% Assign some helpful precomp items
precomp(xDim).Tdif = Tdif;
for i = 1:xDim
    precomp(i).Tdif = Tdif;
    precomp(i).Tall   = Tall;
    precomp(i).params = params;
end
% Find unique numbers of trial lengths
Tu = unique(Tall);
% Loop once for each state dimension (each GP)
for i = 1:xDim
    for j = 1:length(Tu)
        T = Tu(j);
        precomp(i).Tu(j).nList = find(Tall == T);
        precomp(i).Tu(j).T = T;
        precomp(i).Tu(j).numTrials = length(precomp(i).Tu(j).nList);
        precomp(i).Tu(j).Xmoment  = zeros(T*numGroups);
    end
end

% Fill out Xmoment
% Loop once for each state dimension (each GP)
for i = 1:xDim
    % Loop once for each trial length (each of Tu)
    for j = 1:length(Tu)
        % Loop once for each trial (each of nList)
        for n = precomp(i).Tu(j).nList
            xsm_i = seq(n).xsm(i:xDim:xDim*numGroups,:);
            xsm_i = reshape(xsm_i,1,numGroups*Tu(j));
            
            precomp(i).Tu(j).Xmoment = precomp(i).Tu(j).Xmoment +...
                params.X.cov(j).VsmGP(:,:,i) +...
                xsm_i' * xsm_i;
        end
    end
end
