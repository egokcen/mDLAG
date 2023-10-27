% startup.m
%
% Description: This script defines constants used throughout the V1-V2
%              analyses.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Relevant directories

addpath('./data');

%% General constants

% Dataset parameters
binWidth = 20;          % ms

% Processed data files
dataFile = './data/107l003p143_stim1.mat';

% Plotting
units = 'ms';
groupNames = {'V1', 'V2'};
numGroups = length(groupNames);
groupColors = {'#5599FF', '#FF5555'};