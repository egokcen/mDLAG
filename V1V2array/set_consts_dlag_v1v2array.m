% set_consts_dlag_v1v2array.m
%
% Description: This script defines constants for DLAG analyses of the 
%              V1-V2 recordings.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Relevant directories

restoredefaultpath;
addpath('./data');
addpath('./dlag');
addpath(genpath('../DLAG'));

% Access results here
resultDir_dlag = './dlag/results';

%% Fitting parameters

runIdx = 1;
overwriteExisting = true; 
saveData = false;     
method = 'dlag';
numFolds = 0;
fitAll = true;
rGroups = [1 2]; 
xDim_across = 3;       % Selected via cross-validation
xDim_within = [17 1];  % Selected via cross-validation
segLength = 32;       
init_method = 'static'; 
learnDelays = true;   
maxIters = 5e4;         
freqLL = 10;          
freqParam = 100;      
minVarFrac = 0.01;        
randomSeed = 0;
