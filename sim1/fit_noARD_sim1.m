% fit_noARD_sim1.m
%
% Description: This script fits mDLAG models (no ARD) for Simulation 1.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_noARD_sim1;

%% Fit model

% Load synthetic data and ground truth
currDataFile = sprintf('%s/%s', dataDir, dataFile);
ws = load(currDataFile);

% Fit DLAG models
baseDir = sprintf('%s/%s', resultDir, noARD_FitDir);
fit_dlag(runIdx, ws.seqTrue, ...
         'baseDir', baseDir, ...
         'binWidth', binWidth, ...
         'numFolds', numFolds, ...
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
      