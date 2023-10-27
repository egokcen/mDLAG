% set_consts_gfa_sim2.m
%
% Description: This script defines constants for Simulation 2, specific to 
%              group factor analysis (GFA) experiments.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 

%% Setup directories/path

restoredefaultpath;
addpath(dataDir);
addpath(resultDir);
addpath('./util');
addpath(genpath('../gfa'));

%% GFA fitting parameters

xDim_fit = xDim;         
prior_val = 1e-12;    % Set to very small value to set uninformative prior
prior.d.beta = prior_val;
prior.phi.a = prior_val;
prior.phi.b = prior_val;
prior.alpha.a = prior_val;
prior.alpha.b = prior_val;
prior.uv.lbda = 0;
tol = 1e-8;           % Tolerance to determine fitting convergence
maxIters = 3e4;       % Maximum fitting iterations
verbose = true;       % Print fitting progress
randomSeed = 0;       % Set for reproducibility
minVarFrac = 0.0001;
pruneX = true;        % For speed-up, remove latents that become inactive in all groups
saveX = false;        % Set to false to save memory when saving final results
saveCcov = false;     % Set to false to save memory when saving final results
saveFitProgress = true; % Set to true to save lower bound, time each iteration
