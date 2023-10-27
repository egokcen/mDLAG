% startup.m
%
% Description: This script defines constants used throughout Simulation 2:
%              Disentangling concurrent, bidirectional signaling.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
%

%% Relevant directories

dataDir = './data';                   % Access data here
resultDir = './results';              % Access results files here

%% Add directories to matlab path

addpath(dataDir);
addpath(resultDir);
addpath('./util');

%% Synthetic dataset parameters

dataFile = 'data_sim2.mat';    % Name of dataset files

% Dataset size characteristics
N = 1000;                         % Total number of trials
T = 25;                           % Number of samples per sequence
binWidth = 20;                    % Sample period of ground truth data
yDims = [10 10 10];               % Dimensionalities of each observed group
yDim = sum(yDims);
numGroups = length(yDims);        % Total number of groups
xDim = 2;                         % Latent dimensionality
snr = 10*ones(1,numGroups);       % Signal-to-noise ratios of each group
hyperparams.beta = 1;
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;
MAG = 1e6;
hyperparams.a_alpha = MAG.*[1   1;
                            1   1;
                            1   1];
hyperparams.b_alpha = MAG.*ones(numGroups,xDim);
sigDimsTrue = (1./hyperparams.a_alpha) > 0;
units = 'ms';

% Gaussian process (GP) parameters
tau = [50 50];              % GP timescales
eps = 1e-5.*ones(1,xDim);   % GP noise variances
D = [0  0   ;               % Latent delay matrix
     15 -15 ;
     30 -30 ];

%% General constants

groupNames = {'A', 'B', 'C'};
groupColors = {'#5599FF', '#FF5555', '#FF8F00'};