% performance_summary.m
%
% Description: Compare performance of mDLAG to alternative methods.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% mDLAG-0 vs GFA

load('./results/pred_mdlag0.mat');
load('./results/pred_gfa.mat');

figure;
med_height = 8;

subplot(2,2,1);
hold on;
line([0 0.1], [0 0.1], 'color', 'k', 'linestyle', '--');
scatter(R2gfa(:), R2mdlag0(:), 'k.');
axis square;
xlabel('Test R^2, GFA');
ylabel('Test R^2, mDLAG-0');

subplot(2,2,2);
hold on;
line([0 0.15], [0 0.15], 'color', 'k', 'linestyle', '--');
scatter(MSEgfa(:), MSEmdlag0(:), 'k.');
axis square;
xlabel('Test MSE, GFA');
ylabel('Test MSE, mDLAG-0');

subplot(2,2,3);
hold on;
normalization = 'count';
edges = -0.05:0.01:0.05;
line([0 0], [0 med_height], 'color', 'k', 'linestyle', '--');
histogram(R2mdlag0(:) - R2gfa(:), edges,...
      'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', 0.2, 'LineWidth', 2.0, 'Normalization', normalization);
% Plot median
scatter(median(R2mdlag0(:) - R2gfa(:)), med_height, 'vk', 'filled');
xlabel('mDLAG-0 R^2 - GFA R^2');
ylabel(normalization);
axis square;

subplot(2,2,4);
hold on;
normalization = 'count';
edges = -0.008:0.002:0.008;
line([0 0], [0 med_height], 'color', 'k', 'linestyle', '--');
histogram(MSEmdlag0(:) - MSEgfa(:), edges, ...
      'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', 0.2, 'LineWidth', 2.0, 'Normalization', normalization);
% Plot median
scatter(median(MSEmdlag0(:) - MSEgfa(:)), med_height, 'vk', 'filled');
xlabel('mDLAG-0 MSE - GFA MSE');
ylabel(normalization);
axis square;

% Significance tests
[p_R2_mdlag0gfa,~] = signtest(R2mdlag0(:), R2gfa(:), 'alpha', 0.05, 'tail', 'right')

[p_MSE_mdlag0gfa,~] = signtest(MSEmdlag0(:), MSEgfa(:), 'alpha', 0.05, 'tail', 'left')

%% mDLAG vs mDLAG-0

load('./results/pred_mdlag.mat');
load('./results/pred_mdlag0.mat');

figure;
med_height = 8;

subplot(2,2,1);
hold on;
line([0 0.1], [0 0.1], 'color', 'k', 'linestyle', '--');
scatter(R2mdlag0(:), R2mdlag(:), 'k.');
axis square;
xlabel('Test R^2, mDLAG-0');
ylabel('Test R^2, mDLAG');

subplot(2,2,2);
hold on;
line([0 0.15], [0 0.15], 'color', 'k', 'linestyle', '--');
scatter(MSEmdlag0(:), MSEmdlag(:), 'k.')
axis square;
xlabel('Test MSE, mDLAG-0');
ylabel('Test MSE, mDLAG');

subplot(2,2,3);
hold on;
normalization = 'count';
edges = -0.008:0.002:0.008;
line([0 0], [0 med_height], 'color', 'k', 'linestyle', '--');
histogram(R2mdlag(:) - R2mdlag0(:), edges,...
      'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', 0.2, 'LineWidth', 2.0, 'Normalization', normalization);
% Plot median
scatter(median(R2mdlag(:) - R2mdlag0(:)), med_height, 'vk', 'filled');
xlabel('mDLAG R^2 - mDLAG-0 R^2');
ylabel(normalization);
axis square;

subplot(2,2,4);
hold on;
normalization = 'count';
edges = -0.002:0.0004:0.002;
line([0 0], [0 med_height], 'color', 'k', 'linestyle', '--');
histogram(MSEmdlag(:) - MSEmdlag0(:), edges, ...
      'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', 0.2, 'LineWidth', 2.0, 'Normalization', normalization);
% Plot median
scatter(median(MSEmdlag(:) - MSEmdlag0(:)), med_height, 'vk', 'filled');
xlabel('mDLAG MSE - mDLAG-0 MSE');
ylabel(normalization);
axis square;

% Significance tests
[p_R2_mdlag0,~] = signtest(R2mdlag(:), R2mdlag0(:), 'alpha', 0.05, 'tail', 'right')

[p_MSE_mdlag0,~] = signtest(MSEmdlag(:), MSEmdlag0(:), 'alpha', 0.05, 'tail', 'left')

%% mDLAG vs GFA

load('./results/pred_mdlag.mat');
load('./results/pred_gfa.mat');

figure;
med_height = 8;

subplot(2,2,1);
hold on;
line([0 0.1], [0 0.1], 'color', 'k', 'linestyle', '--');
scatter(R2gfa(:), R2mdlag(:), 'k.');
axis square;
xlabel('Test R^2, GFA');
ylabel('Test R^2, mDLAG');

subplot(2,2,2);
hold on;
line([0 0.15], [0 0.15], 'color', 'k', 'linestyle', '--');
scatter(MSEgfa(:), MSEmdlag(:), 'k.');
axis square;
xlabel('Test MSE, GFA');
ylabel('Test MSE, mDLAG');

subplot(2,2,3);
hold on;
normalization = 'count';
edges = -0.05:0.01:0.05;
line([0 0], [0 med_height], 'color', 'k', 'linestyle', '--');
histogram(R2mdlag(:) - R2gfa(:), edges,...
      'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', 0.2, 'LineWidth', 2.0, 'Normalization', normalization);
% Plot median
scatter(median(R2mdlag(:) - R2gfa(:)), med_height, 'vk', 'filled');
xlabel('mDLAG R^2 - GFA R^2');
ylabel(normalization);
axis square;

subplot(2,2,4);
hold on;
normalization = 'count';
edges = -0.008:0.002:0.008;
line([0 0], [0 med_height], 'color', 'k', 'linestyle', '--');
histogram(MSEmdlag(:) - MSEgfa(:), edges, ...
      'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', 0.2, 'LineWidth', 2.0, 'Normalization', normalization);
% Plot median
scatter(median(MSEmdlag(:) - MSEgfa(:)), med_height, 'vk', 'filled');
xlabel('mDLAG MSE - GFA MSE');
ylabel(normalization);
axis square;

% Significance tests
[p_R2_gfa,~] = signtest(R2mdlag(:), R2gfa(:), 'alpha', 0.05, 'tail', 'right')

[p_MSE_gfa,~] = signtest(MSEmdlag(:), MSEgfa(:), 'alpha', 0.05, 'tail', 'left')
