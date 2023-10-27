% mdlag_vs_dlag_v1v2array.m
%
% Description: Compare mDLAG and DLAG performance on the V1-V2 recordings.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% Compare mDLAG and DLAG

load('./mdlag/results/testpred_mdlag.mat');
load('./dlag/results/testpred_dlag.mat');

figure;
subplot(2,2,1);
hold on;
line([0 0.08], [0 0.08], 'color', 'k', 'linestyle', '--');
scatter(R2dlag(:), R2mdlag(:), 'k.');
axis square;
xlabel('Test R^2, DLAG');
ylabel('Test R^2, mDLAG');

subplot(2,2,2);
hold on;
line([0 1], [0 1], 'color', 'k', 'linestyle', '--');
scatter(MSEdlag(:), MSEmdlag(:), 'k.');
axis square;
xlabel('Test MSE, DLAG');
ylabel('Test MSE, mDLAG');

subplot(2,2,3);
hold on;
med_height = 30;
normalization = 'count';
edges = -0.02:0.005:0.02;
line([0 0], [0 med_height], 'color', 'k', 'linestyle', '--');
histogram(R2mdlag(:) - R2dlag(:), edges, ...
      'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', 0.2, 'LineWidth', 2.0, 'Normalization', normalization);
% Plot median
scatter(median(R2mdlag(:) - R2dlag(:)), med_height, 'vk', 'filled');
xlabel('mDLAG R^2 - DLAG R^2');
ylabel(normalization);
axis square;

subplot(2,2,4);
hold on;
med_height = 30;
normalization = 'count';
line([0 0], [0 med_height], 'color', 'k', 'linestyle', '--');
histogram(MSEmdlag(:) - MSEdlag(:), ...
      'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', 0.2, 'LineWidth', 2.0, 'Normalization', normalization);
% Plot median
scatter(median(MSEmdlag(:) - MSEdlag(:)), med_height, 'vk', 'filled');
xlabel('mDLAG MSE - DLAG MSE');
ylabel(normalization);
axis square;

% Significance tests
[p_R2_dlag,~] = signtest(R2mdlag(:), R2dlag(:), 'alpha', 0.05, 'tail', 'right')

[p_MSE_dlag,~] = signtest(MSEmdlag(:), MSEdlag(:), 'alpha', 0.05, 'tail', 'left')