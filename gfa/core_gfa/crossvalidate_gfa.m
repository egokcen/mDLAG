function res = crossvalidate_gfa(Ys, xDim, varargin)
%
% res = crossvalidate_gfa(Ys, xDim, ...)
%
% Description: Cross-validation to determine optimal rank of alpha
%              factorization for group factor analysis.
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
%     numFolds -- int; number of cross-validation folds (default: 4)
%     RList    -- int array; ranks to compare 
%                 (default: 0:min([numGroups xDim])
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
%     tol         -- float; stopping criterion for EM (default: 1e-8)
%     maxIters    -- int; maximum number of EM iterations (default: 1e8)
%     verbose     -- boolean; set true to show CV details (default: false)
%     numWorkers  -- int; Number of cores to use, for parallelization.
%                    (default: 0)
%     randomSeed  -- int or string; seed the random number generator, 
%                    for reproducible initializations (default: 'shuffle')
%     fitAll      -- logical; set to false to avoid fitting a model to all
%                    train data (only relevant if numFolds > 0) 
%                    (default: true)
%     saveX       -- logical; set true to save posterior estimates of
%                    latent variables X. For large datasets, X.mean may be
%                    very large. (default: false)
%     minVarFrac -- float; fraction of overall data variance for each 
%                   observed dimension to set as the private variance 
%                   floor. (default: 0.001)
%
% Outputs:
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
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     02 Jul 2022 -- Initial full revision.
%     09 Jul 2022 -- Added option to save posterior estimates of X.
%     11 Jul 2022 -- Expanded parallelization over both R and folds.
%     12 Jul 2022 -- Added option to remove mean parameters from low-rank
%                    alpha model.

numGroups = length(Ys);
numFolds = 4;
RList = 0:min([numGroups xDim]);
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
randomSeed = 'shuffle';
fitAll = true;
saveX = 'false';
minVarFrac = 0.001;
extraOpts = assignopts(who, varargin);

if ismember(0,RList) && ~include_mu
    fprintf('Error: If R = 0, then mu must be included in the low-rank alpha model.\n');
    res = [];
    return
end

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
rng(randomSeed);
val_indices = [];
if numFolds > 0
    val_indices = crossvalind('Kfold', N, numFolds);
end

% Preallocate output structure
s = struct('R', []);
if numFolds > 0
    s.R2 = [];
    s.R2_sem = [];
    s.MSE = [];
    s.MSE_sem = [];
end
if fitAll
    s.model = [];
end
res = repmat(s,1,length(RList));

% Linearize nested loop, for parallelization
numModels = length(RList)*length(cvFoldList);
% Output structures -- a bit of memory waste, but smooths out parfor usage
R2List = nan(1,numModels);
MSEList = nan(1,numModels);
models = cell(1,numModels);

if numWorkers > 0
    StartParPool(numWorkers);
end

parfor (modelIdx = 1:numModels, numWorkers)
    [j,i] = ind2sub([length(cvFoldList) length(RList)], modelIdx);

    R = RList(i);
    currFold = cvFoldList(j);

    if verbose
        fprintf('Fitting model rank = %d, CV fold %d of %d.\n', R, currFold, numFolds);
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
                 'include_mu', include_mu, ...
                 'randomSeed', randomSeed, ...
                 'saveX', saveX, ...
                 'minVarFrac', minVarFrac);

    if currFold == 0
        % If using all of the data, save parameters
        models{modelIdx} = out;
    else
        % Compute cross-validated leave-group-out performance
        [R2, MSE] = pred_gfa(Yval, out.C, out.d, out.phi);
        R2List(modelIdx) = R2;
        MSEList(modelIdx) = MSE;
    end
    
end

% Save relevant outputs
models = reshape(models,[length(cvFoldList) length(RList)]);
R2List = reshape(R2List,[length(cvFoldList) length(RList)]);
MSEList = reshape(MSEList,[length(cvFoldList) length(RList)]);
for i = 1:length(RList)
    
    res(i).R = RList(i);
    
    if fitAll
        % Only save models fit to all training data
        res(i).model = models{1,i};
    end
    
    if numFolds > 0
        res(i).R2 = mean(R2List(:,i),'omitnan');
        res(i).R2_sem = std(R2List(:,i),'omitnan') / sqrt(numFolds);
        
        res(i).MSE = mean(MSEList(:,i),'omitnan');
        res(i).MSE_sem = std(MSEList(:,i),'omitnan') / sqrt(numFolds);
    end
    
end