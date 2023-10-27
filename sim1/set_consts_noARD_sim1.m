% set_consts_noARD_sim1.m
%
% Description: This script defines constants for Simulation 1, specific to 
%              mDLAG (no ARD) experiments.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Setup directories/path

restoredefaultpath;
addpath(dataDir);
addpath(resultDir);
addpath('./util');
addpath(genpath('../DLAG'));

%% Fitting parameters

runIdx = 1;
overwriteExisting = true; 
saveData = false;
init_method = 'static';       
numFolds = 0;         
rGroups = [1 2];
segLength = Inf;       
learnDelays = true; 
maxIters = 5e4;       
freqLL = 10;          
freqParam = 100;      
minVarFrac = 0.0001;        
randomSeed = 0;
xDim_across = 7;
xDim_within = zeros(1,numGroups);
xDim_cv = 10;