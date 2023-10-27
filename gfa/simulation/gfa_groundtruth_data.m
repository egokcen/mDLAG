% gfa_groundtruth_data.m
%
% Description: This script generates synthetic datasets from the group
%              factor analysis (GFA) model. Ground truth data and 
%              parameters are saved to a file.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 

%% Define GFA ground truth model parameters

rng('shuffle');

% Dataset size characteristics
N = 500;                          % Total number of trials
yDims = [10 10 10];               % Dimensionalities of each observed group
yDim = sum(yDims);                % Joint dimensionality of all groups
M = length(yDims);                % Total number of groups
xDim = 7;                         % Across-group latent dimensionality
snr = 10^(-0.5)*ones(1,M);        % Signal-to-noise ratios of each group

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

% The following are only relevant for low-rank alpha models
R = 'full';      % Rank of alpha model. 'full' leads to a vanilla GFA model
hyperparams.uv_lbda = 1.5;

%% Randomly generate data from a GFA model

data = simdata_gfa(N, yDims, xDim, hyperparams, snr, 'R', R);

%% Save generated groundtruth data and parameters

save('./demo/data/gfa_demo_data_synthetic.mat', 'data', 'hyperparams');