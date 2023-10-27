% data_efficiency_summary.m
%
% Description: Summarize mDLAG test performance as we limit the number of
%              available training trials in an example Neuropixels dataset.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% mDLAG performance as a function of available training trials

NList = [10 25 50 100 150 225];  % No. trials in each training set
stimIdx = 2; % Index into the example dataset
sessIdx = 2; % Index into the example dataset

% Load performances
mdlag_subtrials = load('./results/subtrials_mdlag.mat');
gfa_subtrials = load('./results/pred_gfa.mat');

% mDLAG performance as a function of number of available training trials
figure;
hold on;
plot(NList, mdlag_subtrials.R2,'k.-', 'markersize', 10);
line([0 NList(end)], [gfa_subtrials.R2gfa(stimIdx,sessIdx) ...
    gfa_subtrials.R2gfa(stimIdx,sessIdx)], ...
    'color', 'r', 'linestyle', '--');
ylim([0 0.085]);
axis square;
xlabel('No. training trials');
ylabel('Test R^2');
