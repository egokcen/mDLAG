% mdlag_exampledataset_v1v2array.m
%
% Description: Explore the fit of an mDLAG model to an example V1-V2 
%              dataset.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_v1v2array;

%% Load example dataset and mDLAG results

% Dataset
load(dataFile);
load('./data/traintestsplit.mat');  % Train and test trials

% mDLAG results
res = load(resultFile_mdlag);
res.xDim = res.estParams.xDim;
res.binWidth = binWidth;
seq = dat.seq(test);
estParams = res.estParams;
trackedParams = res.trackedParams;
flags = res.flags;

% Some constants
N = length(seq);

% Check fitting progress
fprintf('\n');
disp(flags)
plotFittingProgress(trackedParams, binWidth, ...
                    'freqLB', freqLB, ...
                    'freqParam', freqParam, ...
                    'units', units);
                
%% Determine dimensionalities of all types

cutoff_sharedvar = 0.02; % Minimum shared variance within a group that a latent must explain
cutoff_snr = 0.001;      % Minimum SNR that a group must have for latents to be signficiant
[dims,sigDims,varExp,dimTypes] = computeDimensionalities(estParams, ...
                                                         cutoff_sharedvar, ...
                                                         cutoff_snr);
                                                     
% Categorize latents
acrossIdxs = find(ismember(sigDims',[1 1],'rows'));
within1Idxs = find(ismember(sigDims',[1 0],'rows'));
within2Idxs = find(ismember(sigDims',[0 1],'rows'));
nsigIdxs = find(ismember(sigDims',[0 0],'rows'));

%% Visualize shared variance explained by each significant latent

% Alpha parameters
alpha_inv_est = 1./estParams.alpha.mean;
% Normalize by the shared variance in each area
alpha_inv_rel_est = alpha_inv_est ./ sum(alpha_inv_est,2);

% Order latents by type and shared variance explained
alpha_across = alpha_inv_rel_est(:,acrossIdxs);
alpha_within1 = alpha_inv_rel_est(:,within1Idxs);
alpha_within2 = alpha_inv_rel_est(:,within2Idxs);
[~, sortA] = sort(alpha_across(2,:), 'descend');
[~, sort1] = sort(alpha_within1(1,:), 'descend');
[~, sort2] = sort(alpha_within2(2,:), 'descend');
acrossIdxs = acrossIdxs(sortA);
within1Idxs = within1Idxs(sort1);
within2Idxs = within2Idxs(sort2);
keptIdxs = [acrossIdxs; within1Idxs; within2Idxs]; % Latents to plot
alpha = [alpha_across(:,sortA) alpha_within1(:,sort1) alpha_within2(:,sort2)];
              
% Plot shared variances
figure;
hold on;
b = bar(alpha','barwidth',1);
for groupIdx = 1:numGroups
    b(groupIdx).FaceColor = groupColors{groupIdx};
end
xlabel('Latent variable');
ylabel('Frac. shared var. exp.');
ylim([0 0.4]);

%% Visualize latent time courses

xspec = 'xve'; % 'xve' gives latent time courses scales by shared variance
               % 'xsm' gives latent time courses with normalized variances

[seqEst,~,~] = inferX(seq, estParams);
seqEst = getSubsetXDims_seq(seqEst,estParams.xDim,numGroups,keptIdxs,'datafield','xsm');
estParamsSub = getSubsetXDims_params(estParams,keptIdxs);
[seqEst, sortParams] = scaleByVarExp(seqEst, estParamsSub, alpha_inv_rel_est(:,keptIdxs), ...
                                     'sortDims', false, ...
                                     'sortGroup', 2, ...
                                     'numDim', 10, ...
                                     'indatafield', 'xsm', ...
                                     'outdatafield', 'xve');
trialIdxs = 51:60;
plotDimsVsTime_mdlag(seqEst(trialIdxs), xspec, sortParams, binWidth, ...
                     'nPlotMax', 10, ...
                     'plotSingle', true, ...
                     'plotMean', false, ...
                     'units', units, ...
                     'trialGroups', {}, ...
                     'trialColors', {});
                 
%% Visualize GP parameters

plotGPparams_mdlag(estParamsSub,binWidth,'sigDims',sigDims(:,keptIdxs),'units',units);
                 
%% Explore leave-group-out prediction performance and SNRs

% Leave-group-out predictive performance
[R2, MSE] = pred_mdlag(seq, estParams);
fprintf('Leave-group-out R^2:\n    %1.4f\n', R2);

% Signal-to-noise ratio of each group, according to estimated mDLAG model
% parameters
snr = computeSNR(estParams.C, estParams.phi);
fprintf('Signal-to-noise ratio of each group:\n    %1.4f  %1.4f\n', snr);
                        