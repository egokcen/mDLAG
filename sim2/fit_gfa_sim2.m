% fit_gfa_sim2.m
%
% Description: This script fits GFA models to data for Simulation 2.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
%

%% Set up script parameters

set_consts_gfa_sim2;

%% Fit GFA model

% Load synthetic data and ground truth
currDataFile = sprintf('%s/%s', dataDir, dataFile);
ws = load(currDataFile);
Ys = seq2cell2D(ws.seqTrue,yDims);
                    
% Fit model
out = em_gfa(Ys, xDim_fit, ...
             'prior', prior, ...
             'tol', tol, ...
             'maxIters', maxIters, ...
             'randomSeed', randomSeed, ...
             'verbose', verbose, ...
             'R', 'full', ...
             'minVarFrac', minVarFrac, ...
             'pruneX', pruneX, ...
             'saveX', saveX, ...
             'saveCcov', saveCcov, ...
             'saveFitProgress', saveFitProgress);

% Save results 
currSaveDir = sprintf('%s', resultDir);
if isfolder(currSaveDir)
    fprintf('Using existing directory %s...\n', currSaveDir);
else
    fprintf('Making directory %s...\n', currSaveDir);
    mkdir(currSaveDir);
end
currSaveFile = sprintf('%s/gfa_results_sim2.mat', currSaveDir);
save(currSaveFile, 'out');
      