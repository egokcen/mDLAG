function res = crossval_pred_gfa(Ys, xDim, varargin)
%
% res = crossval_pred_gfa(Ys, xDim, ...)
%
% Description: Compute cross-validated leave-group-out performance and
%              pairwise prediction performance.
% 
% Arguments:
%
%     Required:
%
%     Ys          -- (1 x numGroups) cell array; list of data matrices 
%                    {(y1Dim x N), (y2Dim x N), ...}
%     xDim        -- int; latent dimensionality
%   
%     Optional:
%
%     R        -- string or int; 'full' for full-rank alpha model, a
%                 positive integer giving the rank of alpha otherwise.
%     numFolds -- int; number of cross-validation folds (default: 4)
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
%     tol         -- float; stopping criterion for EM (default: 1e-8)
%     maxIters    -- int; maximum number of EM iterations (default: 1e8)
%     verbose     -- boolean; set true to show CV details (default: false)
%     numWorkers  -- int; Number of cores to use, for parallelization.
%                    (default: 0)
%     randomSeed  -- int or string; seed the random number generator, 
%                    for reproducible initializations (default: [])
%     fitAll      -- logical; set to false to avoid fitting a model to all
%                    train data (only relevant if numFolds > 0) 
%                    (default: true)
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
%                        iteration time each EM iteration. (default: false)
%
% Outputs:
%     res -- structure with the following fields:
%              R2      -- average cross-validated R^2
%              R2_sem  -- standard error of R2 across CV folds
%              MSE     -- average cross-validated mean-squared error
%              MSE_sem -- standard error of MSE across CV folds
%              model   -- GFA model estimated using all data
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     08 Nov 2022 -- Initial full revision.

numGroups = length(Ys);
numFolds = 4;
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
numWorkers = 0;
randomSeed = [];
fitAll = true;
minVarFrac = 0.001;
pruneX = true;
saveX = false;
saveCcov = false;
saveFitProgress = false;
extraOpts = assignopts(who, varargin);

yDims = zeros(1,numGroups);
for groupIdx = 1:numGroups
    yDims(groupIdx) = size(Ys,1);
end
N = size(Ys{1},2);

if numFolds > 0 && ~fitAll
    cvFoldList = 1:numFolds; % Cross-validation folds, not including training on all data
else
    cvFoldList = 0:numFolds; % Cross-validation folds, including training on all data
end

% Set cross-validation folds (crossvalind randomly generates indices)
if ~isempty(randomSeed)
    rng(randomSeed);
end
val_indices = [];
if numFolds > 0
    val_indices = crossvalind('Kfold', N, numFolds);
end

% Preallocate output structure
res = [];
if numFolds > 0
    res.R2 = [];
    res.R2_sem = [];
    res.MSE = [];
    res.MSE_sem = [];
end
if fitAll
    res.model = [];
end

% Output structures -- a bit of memory waste, but smooths out parfor usage
numModels = length(cvFoldList);
R2_lgo_agg_list = nan(numModels,1);
R2_lgo_indiv_list = nan(numModels,numGroups);
R2_pair_list = nan(numModels,numGroups,numGroups);
R2_diff_list = nan(numModels,numGroups,numGroups);
MSE_lgo_agg_list = nan(numModels,1);
MSE_lgo_indiv_list = nan(numModels,numGroups);
MSE_pair_list = nan(numModels,numGroups,numGroups);
MSE_diff_list = nan(numModels,numGroups,numGroups);
models = cell(1,numModels);

if numWorkers > 0
    StartParPool(numWorkers);
end

for modelIdx = 1:numModels %parfor (modelIdx = 1:numModels, numWorkers)
    
    currFold = cvFoldList(modelIdx);

    if verbose
        fprintf('CV fold %d of %d.\n', currFold, numFolds);
    end

    % Split the data into training and validation sets
    val = false(1, N);
    if currFold > 0
        % Do cross-validation, don't just train on all of the
        % data.
        val = (val_indices == currFold);
    end
    train = ~val;
    Ytrain = cell(1, numGroups);
    Yval = cell(1, numGroups);
    for groupIdx = 1:numGroups
        Ytrain{groupIdx} = Ys{groupIdx}(:,train);
        Yval{groupIdx} = Ys{groupIdx}(:,val);
    end

    % Fit GFA model
    out = em_gfa(Ytrain, xDim, ...
                 'prior', prior, ...
                 'tol', tol, ...
                 'maxIters', maxIters, ...
                 'verbose', false, ...
                 'R', R, ...
                 'pruneX', pruneX, ...
                 'saveX', saveX, ...
                 'saveCcov', saveCcov, ...
                 'saveFitProgress', saveFitProgress);

    if currFold == 0
        % If using all of the data, save parameters
        models{modelIdx} = out;
    else
        % Compute cross-validated leave-group-out performance
        [R2lgo, MSElgo] = pred_gfa(Yval, out.C, out.d, out.phi);
        [R2pair, MSEpair] = pred_gfa_pairwise(Yval, out.C, out.d, out.phi);
        R2diff = R2lgo.indiv - R2pair;
        MSEdiff = MSElgo.indiv - MSEpair;
        
        % Collect outputs
        % R2
        R2_lgo_agg_list(modelIdx) = R2lgo.agg;
        R2_lgo_indiv_list(modelIdx,:) = R2lgo.indiv;
        R2_pair_list(modelIdx,:,:) = R2pair;
        R2_diff_list(modelIdx,:,:) = R2diff;
        % MSE
        MSE_lgo_agg_list(modelIdx) = MSElgo.agg;
        MSE_lgo_indiv_list(modelIdx,:) = MSElgo.indiv;
        MSE_pair_list(modelIdx,:,:) = MSEpair;
        MSE_diff_list(modelIdx,:,:) = MSEdiff;
    end
    
end

% Save relevant outputs
if fitAll
    % Only save models fit to all training data
    res.model = models{1};
end

if numFolds > 0
    
    % R2
    res.R2.lgo.agg = mean(R2_lgo_agg_list, 'omitnan');
    res.R2_sem.lgo.agg = std(R2_lgo_agg_list, 0, 'omitnan') / sqrt(numFolds);
    
    res.R2.lgo.indiv = mean(R2_lgo_indiv_list, 1,'omitnan');
    res.R2_sem.lgo.indiv = std(R2_lgo_indiv_list, 0, 1, 'omitnan') / sqrt(numFolds);
    
    res.R2.pair = squeeze(mean(R2_pair_list, 1, 'omitnan'));
    res.R2_sem.pair = squeeze(std(R2_pair_list, 0, 1, 'omitnan')) / sqrt(numFolds);

    res.R2.diff = squeeze(mean(R2_diff_list, 1, 'omitnan'));
    res.R2_sem.diff = squeeze(std(R2_diff_list, 0, 1, 'omitnan')) / sqrt(numFolds);
    
    % MSE
    res.MSE.lgo.agg = mean(MSE_lgo_agg_list, 'omitnan');
    res.MSE_sem.lgo.agg = std(MSE_lgo_agg_list, 0,'omitnan') / sqrt(numFolds);
    
    res.MSE.lgo.indiv = mean(MSE_lgo_indiv_list, 1, 'omitnan');
    res.MSE_sem.lgo.indiv = std(MSE_lgo_indiv_list, 0, 1, 'omitnan') / sqrt(numFolds);
    
    res.MSE.pair = squeeze(mean(MSE_pair_list, 1, 'omitnan'));
    res.MSE_sem.pair = squeeze(std(MSE_pair_list, 0, 1, 'omitnan')) / sqrt(numFolds);

    res.MSE.diff = squeeze(mean(MSE_diff_list, 1, 'omitnan'));
    res.MSE_sem.diff = squeeze(std(MSE_diff_list, 0, 1, 'omitnan')) / sqrt(numFolds);
end
    