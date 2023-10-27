% ======================================================
% GFA DEMO: Run this script from the main gfa directory
% ======================================================
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% =====================
% 0) Load demo data
% ======================

% Synthetic data generated from a GFA model.
% See simulation/gfa_groundtruth_data.m for a demonstration of how to
% generate simulated data from a GFA model.
dat_file = './demo/data/gfa_demo_data_synthetic';
fprintf('Reading from %s \n',dat_file);
load(dat_file);
numGroups = length(data.Ys);

% Names for each group in the demo data file
groupNames = {'A', 'B', 'C'};
groupColors = {'#5599FF', '#FF5555', '#FF8F00'};

%% =====================
% 1) Fit GFA model
% ======================

xDim = size(data.X,1);
xDim_fit = 10;         % Set to value larger than you think latent dimensionality will be
prior_val = 1e-12;     % Set to very small value to set uninformative prior
prior.d.beta = prior_val;
prior.phi.a = prior_val;
prior.phi.b = prior_val;
prior.alpha.a = prior_val;
prior.alpha.b = prior_val;
prior.uv.lbda = 0;
tol = 1e-8;           % Tolerance to determine fitting convergence
maxIters = 20000;     % Maximum fitting iterations
verbose = true;       % Print fitting progress
randomSeed = 0;       % Set for reproducibility
pruneX = true;        % For speed-up, remove latents that become inactive in all groups
saveX = false;        % Set to false to save memory when saving final results
saveCcov = false;     % Set to false to save memory when saving final results
saveFitProgress = true; % Set to true to save lower bound, time each iteration

out = em_gfa(data.Ys, xDim_fit, ...
             'prior', prior, ...
             'tol', tol, ...
             'maxIters', maxIters, ...
             'randomSeed', randomSeed, ...
             'verbose', verbose, ...
             'pruneX', pruneX, ...
             'saveX', saveX, ...
             'saveCcov', saveCcov, ...
             'saveFitProgress', saveFitProgress);

%% ========================
% 2) Check fitting results
% =========================

% Display flags indicating fitting procedure status
% flags.Convergence -- Indicates that fitting converged according to 'tol'.
%                      Else 'maxIters' was reached before convergence.
% flags.DecreasingLowerBound -- Indicates that the lower bound (objective
%                               function) decreased at some point during
%                               fitting, which suggests there's something
%                               wrong with the fitting or data.
% flags.PrivateVarianceFloor -- Indicates that the variance floor was used 
%                               on one or more observed dimensions. See 
%                               em_gfa.m header for more info.
% flags.xDimsRemoved         -- Number of latent dimensions removed (if 
%                               pruneX is true) due to low variance.
fprintf('\n');
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

%% ==========================================
% 3) Visualize recovery of select parameters
% ===========================================

% Loadings matrices

% Ground truth
Ctrue = vertcat(data.Cs{:});
hinton(Ctrue);

% Estimate
% NOTE: In general, the columns of Cest are unordered, and will not match
%       the order in Ctrue. Here, we can reorder and flip the sign of each
%       dimension to facilitate comparison with the ground truth.
reorder = [4  2  6  7  5  1  3];
rescale = [1 -1  1 -1 -1  1 -1];
Cest = vertcat(out.C.means{:});
hinton(Cest(:,reorder).*rescale);

% Alpha parameters
% The following plots visualize the shared variance explained by each
% latent variable in each area.

% Ground truth
alpha_inv_true = 1./data.alphas;
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
%       the order in alpha_true. Here, we can reorder the latents to 
%       facilitate comparison with the ground truth.
alpha_inv_est = 1./out.alpha.mean;
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

%% ========================================================================
% 4) Explore leave-group-out prediction performance and SNRs on train data
% =========================================================================

% Leave-group-out predictive performance
[R2, MSE] = pred_gfa(data.Ys, out.C, out.d, out.phi);
fprintf('Leave-group-out R^2:\n    %1.4f\n', R2.agg);

% Signal-to-noise ratio of each group, according to estimated GFA model
% parameters
snr = computeSNR(out.C, out.phi);
fprintf('Signal-to-noise ratio of each group:\n    %1.4f  %1.4f  %1.4f\n', snr);

%% ===========================================================
% 5) Determine and then visualize the dimensionalities of all types
% ==================================================================

cutoff_sharedvar = 0.02; % Minimum shared variance within a group that a latent must explain
cutoff_snr = 0.001;      % Minimum SNR that a group must have for latents to be significant
[dims,sigDims,varExp,dimTypes] = computeDimensionalities(out, ...
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

%% =======================================================
% 6) Perform a pairwise analysis of interactions between groups
% ==============================================================

% Compute pairwise dimensionalities and shared variances explained
[pairDims,pairVarExp,pairs] = computeDims_pairs(dims,dimTypes,varExp);

% Visualize results
plotDims_pairs(pairDims, pairs, numGroups, 'groupNames', groupNames);
plotVarExp_pairs(pairVarExp, pairs, numGroups, 'groupNames', groupNames);
