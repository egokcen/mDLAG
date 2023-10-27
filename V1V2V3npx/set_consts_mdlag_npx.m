% set_consts_mdlag_npx.m
%
% Description: This script defines constants for mDLAG analyses of the 
%              Neuropixels recordings.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Relevant directories

restoredefaultpath;
addpath('./data');
addpath('./results');
addpath(genpath('../mDLAG'));

%% Fitting parameters

xDim_fit = 30;           % Number of latent variables to be fitted
randomSeed = 0;          % Seed the random number generator for reproducibility
prior_val = 1e-12;       % Set to very small value to set uninformative prior
prior.d.beta = prior_val;
prior.phi.a = prior_val;
prior.phi.b = prior_val;
prior.alpha.a = prior_val;
prior.alpha.b = prior_val;

tol = 1e-8;           % Tolerance to determine fitting convergence
maxIters = 5e4;       % Maximum fitting iterations
freqLB = 10;          % Check for convergence of lower bound every freqLB iterations
freqParam = 100;      % Store intermediate delay and timescale estimates every freqParam iterations
learnDelays = true;   % Toggle whether to learn delay parameters
verbose = false;      % Print fitting progress
minVarFrac = 0.001;   % Private noise variances will not be allowed to go below this value
maxDelayFrac  = 1.0;  % Maximum delay magnitude (unit: fraction of trial length)
maxTauFrac = 1.0;     % Maximum timescale magnitude (unit: fraction of trial length)
pruneX = true;        % For speed-up, remove latents that become inactive in all groups
saveXcov = false;     % Set to false to save memory when saving final results
saveCcov = false;     % Set to false to save memory when saving final results
