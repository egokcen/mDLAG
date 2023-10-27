% results_summary_noARD_sim1.m
%
% Description: This script displays mDLAG (no ARD) results for 
%              Simulation 1.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_noARD_sim1;

%% Inspect cross-validation results

% Retrieve cross-validated results for all models in the results directory
baseDir = sprintf('%s/%s', resultDir, noARD_CVDir);
[cvResults, ~] = getCrossValResults_dlag(runIdx, 'baseDir', baseDir);

% Plot performance among the candidate models
xDims_grid = zeros(xDim_cv,numGroups+1);
xDims_grid(:,1) = (1:xDim_cv)';
plotPerfvsDim_dlag(cvResults);

%% Load fitted model and ground truth data

currDataFile = sprintf('%s/%s', dataDir, dataFile);
load(dataFile);
baseDir = sprintf('%s/%s', resultDir, noARD_FitDir);
res = getModel_dlag(runIdx, xDim_across, xDim_within, ...
                    'baseDir', baseDir);
                
%% Check fitting results

plotFittingProgress(res, ...
                    'freqLL', freqLL, ...
                    'freqParam', freqParam, ...
                    'units', 'ms');

%% Visualize recovery of the loading matrix

% Loadings matrices

% Ground truth
Ctrue = vertcat(trueParams.Cs{:});
hinton(Ctrue);

% Estimate
% NOTE: In general, the columns of Cest are unordered, and will not match
%       the order in Ctrue. Here, we can reorder the latents to 
%       facilitate comparison with the ground truth.
obs_block_idxs = get_block_idxs(yDims);
lat_block_idxs = get_block_idxs(ones(1,numGroups).*xDim);
Cest = [];

for groupIdx = 1:numGroups
    obsBlock = obs_block_idxs{groupIdx};
    latBlock = lat_block_idxs{groupIdx};
    Cest = [Cest; res.estParams.C(obsBlock(1):obsBlock(2),latBlock(1):latBlock(2))];
end
reorder = [ 4  2  5  3  7  1  6];
rescale = [-1  1  1 -1 -1 -1  1];
hinton(Cest(:,reorder).*rescale);

%% Plot estimated latents

[seqEst, ~] = exactInferenceWithLL_dlag(seqTrue, res.estParams);
plotDimsVsTime_dlag(seqEst, 'xsm', res.estParams, res.binWidth, ...
                  'nPlotMax', 1, ...
                  'plotSingle', false, ...
                  'plotMean', true, ...
                  'units', units);

%% Leave-group-out prediction

[R2, MSE] = pred_dlag(seqTrue, res.estParams);
fprintf('Leave-group-out R^2:\n    %1.4f\n', R2);
              
%% Quantify recovery of latent time courses
        
% Infer latent variables using fitted model
[seqEst, ~] = exactInferenceWithLL_dlag(seqTrue, res.estParams);
Xs_gt = seq2cell2D(seqTrue, trueParams.xDim.*ones(1,numGroups), 'datafield', 'xsm');
Xs_est = seq2cell2D(seqEst, res.estParams.xDim_across.*ones(1,numGroups), 'datafield', 'xsm');

Y_gt = blkdiag(trueParams.Cs{:}) * vertcat(Xs_gt{:}) + vertcat(trueParams.ds{:});
Y_est = res.estParams.C * vertcat(Xs_est{:}) + res.estParams.d;
RSS = sum( sum( ( Y_gt - Y_est ).^2, 1 ) );
TSS = sum( sum( ( Y_gt - repmat( mean(Y_gt,2), [1 size(Y_gt,2)] ) ).^2, 1 ) );
R2_noARD = 1 - RSS / TSS;
fprintf('R^2, latent reconstruction:\n    %1.4f\n', R2_noARD);
