% dlag_exampledataset_v1v2array.m
%
% Description: Explore the fit of a DLAG model to an example V1-V2 
%              dataset.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_dlag_v1v2array;

%% Load example dataset

% Dataset
load(dataFile);
load('./data/traintestsplit.mat');  % Train and test trials

% DLAG results
res = getModel_dlag(runIdx, xDim_across, xDim_within, ...
                    'baseDir', resultDir_dlag);
seq = dat.seq(test);

% Some constants
N = length(seq);
xDim_across = res.estParams.xDim_across;
xDim_within = res.estParams.xDim_within;

% Check fitting progress and GP parameters
plotFittingProgress(res, ...
                    'freqLL', freqLL, ...
                    'freqParam', freqParam, ...
                    'units', 'ms');
                
plotGPparams_dlag(res.estParams, res.binWidth, res.rGroups, ...
                  'plotAcross', true, ...
                  'plotWithin', true, ...
                  'units', 'ms');
               
%% Visualize latent time courses

xspec = 'xve'; % 'xve' gives latent time courses scales by shared variance
               % 'xsm' gives latent time courses with normalized variances

[seqEst, ~] = exactInferenceWithLL_dlag(seq, res.estParams);
total = false;
[varexp, domexp] = computeVarExp_dlag(res.estParams, total);
[seqEst, sortParams] = scaleByVarExp(seqEst, ...
                                     res.estParams, ...
                                     varexp.indiv, ...
                                     'sortDims', true, ...
                                     'sortGroup', 2, ...
                                     'numWithin', 3, ...
                                     'numAcross', 3);
trialIdxs = 51:60;
plotDimsVsTime_dlag(seqEst(trialIdxs), xspec, sortParams, res.binWidth, ...
                  'nPlotMax', 10, ...
                  'plotSingle', true, ...
                  'plotMean', false, ...
                  'units', [], ...
                  'trialGroups', {}, ...
                  'trialColors', {});
              
%% Explore leave-group-out prediction performance

% Leave-group-out predictive performance
[R2, MSE] = pred_dlag(seq, res.estParams);
fprintf('Leave-group-out R^2:\n    %1.4f\n', R2);
