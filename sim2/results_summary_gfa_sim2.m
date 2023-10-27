% results_summary_gfa_sim2.m
%
% Description: This script displays GFA results for Simulation 2.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 

%% Set up script parameters

set_consts_gfa_sim2;

%% Load fitted model and ground truth data

currDataFile = sprintf('%s/%s', dataDir, dataFile);
load(dataFile);
resultFile = sprintf('%s/gfa_results_sim2.mat', resultDir);
load(resultFile);

Ys = seq2cell2D(seqTrue,yDims);

%% Check fitting results

disp(out.flags)

% LOWER BOUND
% The lower bound (objective function) should be monotonically increasing.
figure;
subplot(1,2,1);
hold on;
plot(out.lb, 'k-');
xlabel('Iteration');
ylabel('Lower bound');

% RUNTIME
subplot(1,2,2);
hold on;
plot(cumsum(out.iterTime), 'k-');
xlabel('Iteration');
ylabel('Cumulative runtime (s)');

%% Visualize recovery of the loading matrix and ARD parameters

% Loadings matrices

% Ground truth
Ctrue = vertcat(trueParams.Cs{:});
hinton(Ctrue);

% Estimate
% NOTE: In general, the columns of Cest are unordered, and will not match
%       the order in Ctrue.
Cest = vertcat(out.C.means{:});
hinton(Cest);

%% Visualize recovery of latent time courses

xspec = 'xsm'; % 'xve' gives latent time courses scales by shared variance
               % 'xsm' gives latent time courses with normalized variances 
               
% Ground truth               
plotDimsVsTime(seqTrue, xspec, binWidth, ...
               'nPlotMax', 1, ...
               'nCols', xDim, ...
               'plotSingle', false, ...
               'plotMean', true, ...
               'units', units, ...
               'trialGroups', {}, ...
               'trialColors', {});
           
% Estimates
X = inferX_gfa(Ys, out);
seqEst = dat2seq(reshape(X.mean,[out.xDim T N]), 'datafield', xspec);
plotDimsVsTime(seqEst, xspec, binWidth, ...
               'nPlotMax', 10, ...
               'nCols', xDim, ...
               'plotSingle', true, ...
               'plotMean', false, ...
               'units', units, ...'
               'trialGroups', {}, ...
               'trialColors', {});

%% Explore leave-group-out prediction performance and SNRs on train data

% Leave-group-out predictive performance
[R2, MSE] = pred_gfa(Ys, out.C, out.d, out.phi);
fprintf('Leave-group-out R^2:\n    %1.4f\n', R2.agg);

% Signal-to-noise ratio of each group, according to estimated GFA model
% parameters
snr = computeSNR(out.C, out.phi);
fprintf('Signal-to-noise ratio of each group:\n    %1.4f  %1.4f  %1.4f\n', snr);
           
%% Quantify recovery of latent time courses
        
% Infer latent variables using fitted model
X = inferX_gfa(Ys, out);
Xs_gt = seq2cell2D(seqTrue, trueParams.xDim.*ones(1,numGroups), 'datafield', 'xsm');

% Map latents to observed activity and compare against true (noiseless)
% activity
Y_gt = blkdiag(trueParams.Cs{:}) * vertcat(Xs_gt{:}) + vertcat(trueParams.ds{:});
Y_est = vertcat(out.C.means{:}) * X.mean + out.d.mean;
RSS = sum( sum( ( Y_gt - Y_est ).^2, 1 ) );
TSS = sum( sum( ( Y_gt - repmat( mean(Y_gt,2), [1 size(Y_gt,2)] ) ).^2, 1 ) );
R2_GFA = 1 - RSS / TSS;
fprintf('R^2, latent reconstruction:\n    %1.4f\n', R2_GFA);