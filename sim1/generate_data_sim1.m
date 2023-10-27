% generate_data_demo.m
%
% Description: This script generates synthetic datasets from the mDLAG
%              model. Ground truth data and parameters are saved to files. 
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Randomly generate data from mDLAG model

[seqTrue, trueParams] = simdata_mdlag(N, T, binWidth, yDims, xDim, ...
                                      hyperparams, snr, tau, eps, D);

% Save generated data, along with ground truth parameters
currSaveDir = sprintf('%s', dataDir);
if isfolder(currSaveDir)
    fprintf('Using existing directory %s...\n', currSaveDir);
else
    fprintf('Making directory %s...\n', currSaveDir);
    mkdir(currSaveDir);
end
currSaveFile = sprintf('%s/%s', currSaveDir, dataFile);
save(currSaveFile, 'seqTrue', 'trueParams', 'snr', 'binWidth');

%% Inspect generated data

% GP parameters
gp_params = plotGPparams_mdlag(trueParams,binWidth,'sigDims',sigDimsTrue,'units',units);

% Latent timecourses
% Relative shared variance explained by each dimension
alpha_inv = 1./trueParams.alphas;
alpha_inv_rel = alpha_inv ./ sum(alpha_inv,2);
[seqTrue, sortParams] = scaleByVarExp(seqTrue, trueParams, alpha_inv_rel, ...
                                      'sortDims', false, ...
                                      'sortGroup', 1, ...
                                      'numDim', 7, ...
                                      'indatafield', 'xsm', ...
                                      'outdatafield', 'xve');
xspec = 'xve';
plotDimsVsTime_mdlag(seqTrue, xspec, sortParams, binWidth, ...
                     'nPlotMax', 10, ...
                     'plotSingle', true, ...
                     'plotMean', false, ...
                     'units', units, ...
                     'trialGroups', {}, ...
                     'trialColors', {});
