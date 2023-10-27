function Xmoment = computeXmoment(seq,params)
%
% Xmoment = computeXmoment(seq,params)
%
% Description: Helper function to compute the posterior second moments
%              of X for each group.
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
%
%     Xmoment -- (1 x numGroups) cell array; Xmoment{m} is a (xDim x xDim)
%                array containing second moments of X for the mth group.           
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

% Find unique numbers of trial lengths
Tu = unique(Tall);

% Fill out Xmoment
Xmoment = cell(1,numGroups);
for groupIdx = 1:numGroups
    % Get the appropriate latent indices for the current group
    lat_idxs = (1+(groupIdx-1)*xDim):(groupIdx*xDim);
    Xmoment{groupIdx} = zeros(xDim);
    % Loop over each trial length (each of Tu)
    for j = 1:length(Tu)
        % Loop  over each trial (each of nList)
        nList = find(Tall == Tu(j));
        for n = nList
            Xmoment{groupIdx} = Xmoment{groupIdx} + ...
                sum(params.X.cov(j).Vsm(lat_idxs,lat_idxs,:), 3) ...
                    + seq(n).xsm(lat_idxs,:) * seq(n).xsm(lat_idxs,:)';
        end
    end
end
