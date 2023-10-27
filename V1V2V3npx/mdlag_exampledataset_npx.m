% mdlag_exampledataset_npx.m
%
% Description: Explore the fit of an mDLAG model to an example Neuropixels 
%              dataset.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_npx;

%% Load example model and outputs

load('./results/model_mdlag.mat');
load('./results/examplelatents_mdlag.mat');

%% Fitting progress

disp(model.flags)
plotFittingProgress(model.trackedParams, binWidth, ...
                    'freqLB', freqLB, ...
                    'freqParam', freqParam, ...
                    'units', units);
     
%% Dimensionalities

cutoff_sharedvar = 0.02; % Minimum shared variance within a group that a latent must explain
cutoff_snr = 0.001;      % Minimum SNR that a group must have for latents to be signficiant
[dims,sigDims,varExp,dimTypes] = computeDimensionalities(model.estParams, ...
                                                         cutoff_sharedvar, ...
                                                         cutoff_snr);                                        

% Identify global vs unique pairwise interactions
globalIdxs = find(ismember(sigDims',[1 1 1],'rows'));
v1v2Idxs = find(ismember(sigDims',[1 1 0],'rows'));
v1v3Idxs = find(ismember(sigDims',[1 0 1],'rows'));
v2v3Idxs = find(ismember(sigDims',[0 1 1],'rows'));
keptIdxs = [globalIdxs; v1v2Idxs; v1v3Idxs; v2v3Idxs];

%% Alpha parameters
% The following plots visualize the shared variance explained by each
% latent variable in each area.
alpha_inv = 1./model.estParams.alpha.mean;
% Normalize by the shared variance in each area
alpha_inv_rel = alpha_inv ./ sum(alpha_inv,2);

% Order according to interaction type
alpha_global = alpha_inv_rel(:,globalIdxs);
alpha_v1v2 = alpha_inv_rel(:,v1v2Idxs);
alpha_v1v3 = alpha_inv_rel(:,v1v3Idxs);
alpha_v2v3 = alpha_inv_rel(:,v2v3Idxs);
alpha = [alpha_global alpha_v1v2 alpha_v1v3 alpha_v2v3];

% Plot shared variances
figure;
hold on;
b = bar(alpha','barwidth',1);
for groupIdx = 1:numGroups
    b(groupIdx).FaceColor = groupColors{groupIdx};
end
xlabel('Latent variable');
ylabel('Frac. shared var. exp.');
ylim([0 0.6]);

%% Latent time courses
         
seqSub = getSubsetXDims_seq(seqEst,model.estParams.xDim,numGroups,keptIdxs,'datafield','xsm');
estParamsSub = getSubsetXDims_params(model.estParams,keptIdxs);
plotDimsVsTime_mdlag(seqSub, 'xsm', estParamsSub, binWidth, ...
                     'nPlotMax', 10, ...
                     'plotSingle', true, ...
                     'plotMean', false, ...
                     'plotErr', 0, ...
                     'units', units, ...
                     'trialGroups', {}, ...
                     'trialColors', {});