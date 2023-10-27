% startup.m
%
% Description: This script defines constants used throughout the
%              Neuropixels analyses.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Relevant directories

addpath ./data
addpath ./results

%% General constants

% Dataset parameters
binWidth = 20;
numStim = 2;
numSess = 5;

% Plotting
units = 'ms';
groupNames = {'V1', 'V2', 'V3d'};
numGroups = length(groupNames);
groupColors = {'#5599FF', '#FF5555', '#FF8F00'};
