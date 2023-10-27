function bestModel = getBestModelRank_gfa(res, varargin)
%
% bestModel = getBestModelRank_gfa(res, ...)
%
% Description: Determine the best model rank given the cross-validated
%              results in res.
% 
% Arguments:
%
%     Required:
%
%     res -- structure whose i-th entry (corresponding to the i-th latent
%            dimensionality) has fields
%              R       -- rank of factorization
%              R2      -- average cross-validated leave-group-out R^2
%              R2_sem  -- standard error of R2 across CV folds
%              MSE     -- average cross-validated leave-group-out 
%                         mean-squared error
%              MSE_sem -- standard error of MSE across CV folds
%              model   -- GFA model estimated using all data
%
%     Optional:
%     
%     selection   -- string; metric used to select best model ('R2' or
%                    'MSE'; default: 'R2')
%
% Outputs:
%
%     bestModel -- int; index corresponding to the best model in res
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     02 Jul 2022 -- Initial full revision.

selection = 'R2';
extraOpts = assignopts(who, varargin);
    
% Optimal dimensionality according to smallest value whose 
% performance is within one SEM of the best performance
if isequal(selection, 'R2')
    R2 = [res.R2];
    [best_performance, maxIdx] = max(R2);
    best_sem = res(maxIdx).R2_sem;
    % NOTE: max used in this manner selects the first instance that
    %       satisfies the input condition.
    [~, bestModel] = max(R2 >= best_performance - best_sem); 
else
    MSE = [res.MSE];
    [best_performance, minIdx] = min(MSE);
    best_sem = res(minIdx).MSE_sem;
    % NOTE: min used in this manner selects the first instance that
    %       satisfies the input condition.
    [~, bestModel] = min(MSE <= best_performance + best_sem); 
end