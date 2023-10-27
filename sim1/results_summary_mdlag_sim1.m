% results_summary_mdlag_sim1.m
%
% Description: This script displays mDLAG results for Simulation 1.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_sim1;

%% Load fitted model and ground truth data

currDataFile = sprintf('%s/%s', dataDir, dataFile);
load(dataFile);
resultFile = sprintf('%s/mdlag_results_sim1.mat', resultDir);
load(resultFile);

%% Check fitting results

disp(flags)
plotFittingProgress(trackedParams, binWidth, ...
                    'freqLB', freqLB, ...
                    'freqParam', freqParam, ...
                    'units', units);

%% Visualize recovery of the loading matrix and ARD parameters

% Loadings matrices

% Ground truth
Ctrue = vertcat(trueParams.Cs{:});
hinton(Ctrue);

% Estimate
% NOTE: In general, the columns of Cest are unordered, and will not match
%       the order in Ctrue. Here, we can reorder the latents to 
%       facilitate comparison with the ground truth.
Cest = vertcat(estParams.C.means{:});
reorder = [ 3  2  5 6 1  4 7];
rescale = [-1 -1 -1 1 1 -1 1];
hinton(Cest(:,reorder).*rescale);

% Alpha parameters
% The following plots visualize the shared variance explained by each
% latent variable in each area.

% Ground truth
alpha_inv_true = 1./trueParams.alphas;
% Normalize by the shared variance in each area
alpha_inv_rel_true = alpha_inv_true ./ sum(alpha_inv_true,2);

figure;
hold on;
b = bar(alpha_inv_rel_true');
for groupIdx = 1:numGroups
    b(groupIdx).FaceColor = groupColors{groupIdx};
end
xlabel('Latent variable');
ylabel('Frac. shared var. exp.');
ylim([0 0.5]);
title('Ground truth')

% Estimate
% NOTE: In general, the columns of alpha_est are unordered, and will not match
%       the order in alpha_true.
alpha_inv_est = 1./estParams.alpha.mean;
% Normalize by the shared variance in each area
alpha_inv_rel_est = alpha_inv_est ./ sum(alpha_inv_est,2);

figure;
hold on;
b = bar(alpha_inv_rel_est(:,reorder)');
for groupIdx = 1:numGroups
    b(groupIdx).FaceColor = groupColors{groupIdx};
end
xlabel('Latent variable');
ylabel('Frac. shared var. exp.');
ylim([0 0.5]);
title('Estimate')

%% Determine and then visualize the dimensionalities of all types

cutoff_sharedvar = 0.02; % Minimum shared variance within a group that a latent must explain
cutoff_snr = 0.001;      % Minimum SNR that a group must have for latents to be significant
[dims,sigDims,varExp,dimTypes] = computeDimensionalities(estParams, ...
                                                         cutoff_sharedvar, ...
                                                         cutoff_snr);
                                              
% Visualize the number of each type of dimension
plotDimensionalities(dims, dimTypes, ...
                     'groupNames', groupNames, ...
                     'plotZeroDim', false);

% Visualize the shared variance explained by each dimension type in each
% group
plotVarExp(varExp, dimTypes, ...
           'groupNames', groupNames, ...
           'plotZeroDim', false);

%% Visualize recovery of GP parameters

% Ground truth
plotGPparams_mdlag(trueParams,binWidth,'sigDims',sigDimsTrue,'units',units);

% Estimate
plotGPparams_mdlag(estParams,binWidth,'sigDims',sigDims,'units',units);

%% Visualize recovery of latent time courses

xspec = 'xsm'; % 'xve' gives latent time courses scales by shared variance
               % 'xsm' gives latent time courses with normalized variances
                 
% Ground truth
[seqTrue, sortParams] = scaleByVarExp(seqTrue, trueParams, alpha_inv_rel_true, ...
                                     'sortDims', false, ...
                                     'sortGroup', 1, ...
                                     'numDim', 10, ...
                                     'indatafield', 'xsm', ...
                                     'outdatafield', 'xve');
plotDimsVsTime_mdlag(seqTrue, xspec, sortParams, binWidth, ...
                     'nPlotMax', 1, ...
                     'plotSingle', false, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {}, ...
                     'trialColors', {});
                 
% Estimate
[seqEst,~,~] = inferX(seqTrue, estParams);
% NOTE: In general, latent time courses are unordered, and will not match
%       the order in seqTrue.
[seqEst, sortParams] = scaleByVarExp(seqEst, estParams, alpha_inv_rel_est, ...
                                     'sortDims', false, ...
                                     'sortGroup', 1, ...
                                     'numDim', 10, ...
                                     'indatafield', 'xsm', ...
                                     'outdatafield', 'xve');
plotDimsVsTime_mdlag(seqEst, xspec, sortParams, binWidth, ...
                     'nPlotMax', 1, ...
                     'plotSingle', false, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {}, ...
                     'trialColors', {});
                 
% Overlay estimate on ground truth, for an example trial                 
trialIdx = 1;
seqEst(trialIdx).(xspec) = seqEst(trialIdx).(xspec)([reorder reorder+estParams.xDim reorder+estParams.xDim*2],:) .* repmat(rescale, 1, 3)';                
plotDimsVsTime_mdlag([seqTrue(trialIdx); seqEst(trialIdx)], xspec, trueParams, binWidth, ...
                     'nPlotMax', 1, ...
                     'plotSingle', false, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {[1], [2]}, ...
                     'trialColors', {'k', '#D35FBC'});

%% Explore leave-group-out prediction performance and SNRs on train data

% Leave-group-out predictive performance
[R2, MSE] = pred_mdlag(seqTrue, estParams);
fprintf('Leave-group-out R^2:\n    %1.4f\n', R2);

% Signal-to-noise ratio of each group, according to estimated mDLAG model
% parameters
snr = computeSNR(estParams.C, estParams.phi);
fprintf('Signal-to-noise ratio of each group:\n    %1.4f  %1.4f  %1.4f\n', snr);
                 
%% Quantify recovery of latent time courses
        
% Infer latent variables using fitted model
[seqEst,~,~] = inferX(seqTrue, estParams);
Xs_gt = seq2cell2D(seqTrue, trueParams.xDim.*ones(1,numGroups), 'datafield', 'xsm');
Xs_est = seq2cell2D(seqEst, estParams.xDim.*ones(1,numGroups), 'datafield', 'xsm');

Y_gt = blkdiag(trueParams.Cs{:}) * vertcat(Xs_gt{:}) + vertcat(trueParams.ds{:});
Y_est = blkdiag(estParams.C.means{:}) * vertcat(Xs_est{:}) + estParams.d.mean;
RSS = sum( sum( ( Y_gt - Y_est ).^2, 1 ) );
TSS = sum( sum( ( Y_gt - repmat( mean(Y_gt,2), [1 size(Y_gt,2)] ) ).^2, 1 ) );
R2_mDLAG = 1 - RSS / TSS;
fprintf('R^2, latent reconstruction:\n    %1.4f\n', R2_mDLAG);
