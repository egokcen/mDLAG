% results_summary_mdlag_sim2.m
%
% Description: This script displays mDLAG results for Simulation 2.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 

%% Set up script parameters

set_consts_mdlag_sim2;

%% Load fitted model and ground truth data

currDataFile = sprintf('%s/%s', dataDir, dataFile);
load(dataFile);
resultFile = sprintf('%s/mdlag_results_sim2.mat', resultDir);
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
reorder = [ 2  1];
rescale = [ 1 -1];
hinton(Cest(:,reorder).*rescale);

%% Visualize recovery of latent time courses
                 
% Ground truth
plotDimsVsTime_mdlag(seqTrue, 'xsm', trueParams, binWidth, ...
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
plotDimsVsTime_mdlag(seqEst, 'xsm', estParams, binWidth, ...
                     'nPlotMax', 10, ...
                     'plotSingle', true, ...
                     'plotMean', false, ...
                     'units', units, ...
                     'trialGroups', {}, ...
                     'trialColors', {});

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