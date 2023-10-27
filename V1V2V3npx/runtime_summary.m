% runtime_summary.m
%
% Description: Summarize mDLAG runtime across Neuropixels datasets and
%              compare to GFA runtimes.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% Compare mDLAG runtimes to GFA

% Load data characteristics and runtimes
load('./data/population_sizes.mat');
mdlag_runtimes = load('./results/runtime_mdlag.mat');
gfa_runtimes = load('./results/runtime_gfa.mat');

% Average runtime per iteration vs number of neurons in each session
figure;
hold on;
for sessIdx = 1:numSess
    for stimIdx = 1:numStim
        plot(sum(pop_sizes{sessIdx}), mdlag_runtimes.itertimes(stimIdx,sessIdx), 'k.');
        plot(sum(pop_sizes{sessIdx}), gfa_runtimes.itertimes(stimIdx,sessIdx), '.', ...
             'color', [0.5 0.5 0.5]);
    end
end
xlabel('Total no. neurons');
ylabel('Avg. runtime per iter. (s)');
axis square;
