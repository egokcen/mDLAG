% generate_data_sim2.m
%
% Description: This script generates synthetic datasets for Simulation 2:
%              Disentangling concurrent, bidirectional signaling.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_sim2;

%% Generate data

% Generate initial mDLAG model parameters
trueParams = generate_params_mdlag(yDims, xDim, binWidth, hyperparams, snr, ...
                                   tau, eps, D);
                           
% Generate initial latent sequences
seqTrue = generate_latents_mdlag(trueParams, T, N, 'latentfield', 'xsm');
 
% Modify dimensions to have equal strength
for groupIdx = 1:numGroups
    trueParams.Cs{groupIdx}(:,sigDimsTrue(groupIdx,:)) = normc(trueParams.Cs{groupIdx}(:,sigDimsTrue(groupIdx,:)));
    % Enforce the desired signal-to-noise ratios
    varCC = trace(trueParams.Cs{groupIdx} * trueParams.Cs{groupIdx}');
    varNoise_desired = varCC / snr(groupIdx);
    varNoise_current = sum(trueParams.phis{groupIdx}.^(-1));
    trueParams.phis{groupIdx} = trueParams.phis{groupIdx} .* (varNoise_current / varNoise_desired);
end

% Modify latents to have recognizable 'hills' of activity
seqLength = T*binWidth; % Length, in ms
peaks = [150 350]; % ms
timescales = [40 40];
offset = [150 -100];
for n = 1:N
    for j = 1:length(timescales)
        for groupIdx = 1:numGroups
            xIdx = j+(groupIdx-1)*xDim;
            seqTrue(n).xsm(xIdx,:) = 1.*exp(-0.5*(1/timescales(j)).^2*(((1:binWidth:seqLength)-peaks(j))-D(groupIdx,j)).^2) ...
                - 0.*exp(-0.5*(1/timescales(j)).^2*(((1:binWidth:seqLength)-peaks(j))-D(groupIdx,j)-offset(j)).^2); 
        end
    end
end

% Generate observed sequences
seqTrue = generate_obs_mdlag(seqTrue, trueParams, 'latentfield', 'xsm', 'obsfield', 'y');

% Save generated data, along with ground truth parameters
currSaveDir = sprintf('%s', dataDir);
if isfolder(currSaveDir)
    fprintf('Using existing directory %s...\n', currSaveDir);
else
    fprintf('Making directory %s...\n', currSaveDir);
    mkdir(currSaveDir);
end
currSaveFile = sprintf('%s/%s', currSaveDir, dataFile);
save(currSaveFile, 'seqTrue', 'trueParams', 'snr', 'binWidth');

%% Inspect generated data

% GP parameters
gp_params = plotGPparams_mdlag(trueParams,binWidth,'sigDims',sigDimsTrue,'units',units);

xspec = 'xsm';
plotDimsVsTime_mdlag(seqTrue, xspec, trueParams, binWidth, ...
                     'nPlotMax', 10, ...
                     'plotSingle', true, ...
                     'plotMean', false, ...
                     'units', units, ...
                     'trialGroups', {}, ...
                     'trialColors', {});
