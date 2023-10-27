% mdlag_groundtruth_data.m
%
% Description: This script generates synthetic datasets from the mDLAG
%              model. Ground truth data and parameters are saved to a file.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
%

%% Define mDLAG ground truth model parameters

rng('shuffle');

% Dataset size characteristics
N = 100;                          % Number of sequences (trials)
T = 25;                           % Number of samples per sequence
binWidth = 20;                    % Sample period of ground truth data
yDims = [10 10 10];               % Dimensionalities of each observed group
yDim = sum(yDims);                % Joint dimensionality of all groups
M = length(yDims);                % Total number of groups
xDim = 7;                         % Across-group latent dimensionality
snr = 1.0*ones(1,M);              % Signal-to-noise ratios of each group

% a_alpha (M x xDim array) is where the inter-group interaction structure 
% is defined. Row i corresponds to group i. Column j corresponds to latent 
% j. A value of Inf indicates that a latent is NOT present in a group. The
% corresponding loadings will be 0 for that group.
MAG = 20; % Control the variance of alpha parameters (larger = less var.)
hyperparams.a_alpha = MAG.*[1 1   1   Inf 1   Inf Inf;
                            1 1   Inf 1   Inf 1   Inf;
                            1 Inf 1   1   Inf Inf 1  ];
% The remaining hyperparameters are not very important, and can be left
% alone.
hyperparams.b_alpha = MAG.*ones(M,xDim);
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;
hyperparams.beta = 1;

% Gaussian process (GP) parameters
tau = [30 80 50 120 100 40 70];   % GP timescales
eps = 1e-3.*ones(1,xDim);         % GP noise variances
D = [0    0   0  0  0  0  0;      % Latent delay matrix
     15 -30   0  0  0  0  0;
     30   0 -25 40  0  0  0];

%% Randomly generate data from a mDLAG model

[seqTrue, paramsTrue] = simdata_mdlag(N, T, binWidth, yDims, xDim, hyperparams, snr, tau, eps, D);

%% Visualize the ground truth

% GP parameters
sigDimsTrue = (1./hyperparams.a_alpha) > 0;
units = 'ms';
gp_params = plotGPparams_mdlag(paramsTrue,binWidth,'sigDims',sigDimsTrue,'units',units);

% Latent timecourses
% Relative shared variance explained by each dimension
alpha_inv = 1./paramsTrue.alphas;
alpha_inv_rel = alpha_inv ./ sum(alpha_inv,2);
[seqTrue, sortParams] = scaleByVarExp(seqTrue, paramsTrue, alpha_inv_rel, ...
                                      'sortDims', false, ...
                                      'sortGroup', 1, ...
                                      'numDim', 7, ...
                                      'indatafield', 'xsm', ...
                                      'outdatafield', 'xve');
xspec = 'xve';
plotDimsVsTime_mdlag(seqTrue, xspec, sortParams, binWidth, ...
                     'nPlotMax', 10, ...
                     'plotSingle', true, ...
                     'plotMean', false, ...
                     'units', units, ...
                     'trialGroups', {}, ...
                     'trialColors', {});

%% Save generated ground truth data and parameters

save('./demo/data/mdlag_demo_data_synthetic.mat', 'seqTrue', 'paramsTrue', 'binWidth', 'sigDimsTrue');