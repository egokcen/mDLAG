function Kpred = construct_K_pred(params,T,predGroup,varargin)
%
% Kpred = construct_K_pred(params,T,predGroup,...)
%
% Description: Construct portion of GP covariance matrix that can be used
%              to predict predGroup's latent variables from the latent
%              variables of the remaining groups.
%
% Arguments:
%
%     Required:
%
%     params    -- Structure containing mDLAG model parameters.
%                  Contains the following relevant fields:
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
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
%
%     T         -- int; number of timesteps
%     predGroup -- int; index of group to be predicted
%
%     Optional:
%
%     sourceGroups -- (1 x numSourceGroups) array; indices of desired
%                     groups. By default, sourceGroups is the complement
%                     of predGroup, i.e., it includes all groups BUT
%                     predGroup.
%
% Outputs:
%
%     Kpred -- (xDim * T) x (xDim * (numGroups-1) * T) array; predictive GP 
%              covariance matrix            
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     22 Oct 2022 -- Initial full revision.

xDim         = params.xDim;
yDims        = params.yDims;
numGroups    = length(yDims);
sourceGroups = setdiff(1:numGroups,predGroup);
assignopts(who,varargin);

% Construct full GP covariance matrix
K = construct_K_mdlag(params,T);

% Construct a TxT cell array of numGroups*xDim x numGroups*xDim array
% blocks
K = mat2cell(K,repmat(numGroups*xDim,1,T),repmat(numGroups*xDim,1,T));

Kts = cell(T,T); % GP covariance, source-to-target
Ks = cell(T,T);  % GP covariance, source-to-source
for t1 = 1:T
    for t2 = 1:T
        % Convert each array block of K into a numGroups x numGroups cell
        % array
        K{t1,t2} = mat2cell(K{t1,t2},repmat(xDim,1,numGroups),repmat(xDim,1,numGroups));
        % Then, keep only rows corresponding to predGroup, and the columns
        % corresponding to sourceGroups
        Kts{t1,t2} = K{t1,t2}(predGroup,sourceGroups);
        % Convert back to array
        Kts{t1,t2} = cell2mat(Kts{t1,t2});
        
        % Next, keep only rows and columns corresponding to sourceGroups 
        Ks{t1,t2} = K{t1,t2}(sourceGroups,sourceGroups);
        % Convert back to array
        Ks{t1,t2} = cell2mat(Ks{t1,t2});
    end
end

% Convert Kts and Ks back to arrays
Kts = cell2mat(Kts);
Ks = cell2mat(Ks);

% Compute Kpred
Kpred = Kts/Ks;
