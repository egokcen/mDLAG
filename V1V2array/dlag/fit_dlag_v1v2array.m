% fit_dlag_v1v2array.m
%
% Description: Fit a DLAG model to the example V1-V2 dataset.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_dlag_v1v2array;

%% Fit DLAG model

% Load dataset
load(dataFile);
load('./data/traintestsplit.mat');  % Train and test trials
yDims = dat.yDims;

fit_dlag(runIdx, dat.seq(train), ...
         'baseDir', resultDir_dlag, ...
         'method', method, ...
         'binWidth', binWidth, ...
         'numFolds', numFolds, ...
         'fitAll', fitAll, ...
         'xDims_across', xDim_across, ...
         'xDims_within', num2cell(xDim_within), ...
         'yDims', yDims, ...
         'rGroups', rGroups,...
         'segLength', segLength, ...
         'init_method', init_method, ...
         'learnDelays', learnDelays, ...
         'maxIters', maxIters, ...
         'freqLL', freqLL, ...
         'freqParam', freqParam, ...
         'minVarFrac', minVarFrac, ...
         'parallelize', false, ...
         'randomSeed', randomSeed, ...
         'overwriteExisting', overwriteExisting, ...
         'saveData', saveData);