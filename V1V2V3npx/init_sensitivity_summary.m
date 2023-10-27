% init_sensitivity_summary.m
%
% Description: Summarize the sensitivity of mDLAG Neuropixels results to
%              random initialization.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_npx;

%% Collect results and dimensionalities across runs

numRuns = 20;

cutoff_sharedvar = 0.02;
cutoff_snr = 0.001;

% Load results
load('./results/init_sensitivity_mdlag.mat');
numGroups = 3;
groupNames = groupNames(1:numGroups);
numDimTypes = 2^numGroups;
dims = nan(numRuns,numDimTypes);
dimTypes = generateDimTypes(numGroups);
snrs = nan(numGroups,numRuns);
varExp = nan(numRuns,numGroups,numDimTypes);
pairs = nchoosek(1:numGroups,2);
numPairs = size(pairs,1);
pairDims = nan(numRuns,numPairs,3);
pairVarExp = nan(numRuns,numPairs,2);

for runIdx = 1:numRuns

    model = results{runIdx};
    snrEst = computeSNR(model.estParams.C, model.estParams.phi);
    snrs(:,runIdx) = snrEst;

    % All dimension types
    [dims(runIdx,:),sigdims,varExp(runIdx,:,:),~] ...
        = computeDimensionalities(model.estParams, ...
                                  cutoff_sharedvar, ...
                                  cutoff_snr);

    % Pairwise analysis
    [pairDims(runIdx,:,:), pairVarExp(runIdx,:,:), ~] ...
        = computeDims_pairs(squeeze(dims(runIdx,:)), dimTypes, squeeze(varExp(runIdx,:,:)));

    % Collect GP parameters
    if runIdx == 1
        jointParams.yDims = model.estParams.yDims;
        jointParams.xDim = model.estParams.xDim;
        jointParams.D = model.estParams.D;
        jointParams.gamma = model.estParams.gamma;
        jointParams.sigdims = sigdims;
        exampleRun = jointParams;
    else
        jointParams.xDim = jointParams.xDim + model.estParams.xDim;
        jointParams.D = [jointParams.D model.estParams.D];
        jointParams.gamma = [jointParams.gamma model.estParams.gamma];
        jointParams.sigdims = [jointParams.sigdims sigdims];
    end
end

%% Overlay total global on total shared dimensions, separate runs

% Set up axis labels
lbls = cell(1,numPairs);
for pairIdx = 1:numPairs
    % Shared between both groups
    lbls{pairIdx} = sprintf('%s-%s', groupNames{pairs(pairIdx,1)}, groupNames{pairs(pairIdx,2)});
end
% Plot dimensionalities
% Session by session
figure;
for runIdx = 1:numRuns
    subplot(1,numRuns,runIdx);
    hold on;
    for pairIdx = 1:numPairs

        totalShared = pairDims(runIdx,pairIdx,2);
        globalDims = dims(runIdx,end);
        bar(pairIdx, totalShared,  ...
            'facealpha', 0.2, ...
            'facecolor', 'k', ...
            'edgecolor', 'k', ...
            'linewidth', 1.5);
        bar(pairIdx, globalDims,  ...
            'facealpha', 0.5, ...
            'facecolor', 'k', ...
            'edgecolor', 'k', ...
            'linewidth', 1.5);

    end
    ylabel('No. shared dimensions');
    axis([0.5 3.5 0 10]);
    xticklabels(lbls);
    title(sprintf('Run %d', runIdx));
end

%% Plot timescales and delays

% Identify global and unique interactions
noGlobal = jointParams.sigdims;
noGlobal(:,ismember(jointParams.sigdims',ones(1,numGroups),'rows')) = 0;
onlyGlobal = jointParams.sigdims;
onlyGlobal(:,~ismember(jointParams.sigdims',ones(1,numGroups),'rows')) = 0;

% Single out parameters for example run
noGlobal_example = exampleRun.sigdims;
noGlobal_example(:,ismember(exampleRun.sigdims',ones(1,numGroups),'rows')) = 0;
onlyGlobal_example = exampleRun.sigdims;
onlyGlobal_example(:,~ismember(exampleRun.sigdims',ones(1,numGroups),'rows')) = 0;
              
% Convert GP params into units of time
gp_params = getGPparams_mdlag(jointParams, binWidth);
gp_params_example = getGPparams_mdlag(exampleRun, binWidth);

figure;
maxTau = 140;
maxDelay = 35;
colors = generateColors(); % Generate custom plotting colors
pointsize = 5; % Size of scatterplot points
for i = 1:numGroups
    for j = 1:numGroups

        % Consider only different groups
        if i ~= j
            % Take only the delays associated with the current pair of groups
            delays = gp_params.D(j,:) ...
                   - gp_params.D(i,:);
            delays_example = gp_params_example.D(j,:) ...
                           - gp_params_example.D(i,:);

            % Determine which global latents are shared between the current pair of
            % groups
            sig = find(all(onlyGlobal([i j],:),1));
            sig_example = find(all(onlyGlobal_example([i j],:),1));

            subplot(numGroups,numGroups,(i-1)*numGroups+j);
            hold on;
            xlabel(sprintf('Delay from %s to %s (%s)',groupNames{i},groupNames{j},units));
            ylabel(sprintf('GP timescale %s', units));
            xlim([-1.1*maxDelay,1.1*maxDelay]);
            ylim([0,1.1*maxTau]);
            line([0 0], [0 1.1*maxTau], ...
                 'Color', colors.grays{6}, ...
                 'linestyle', '--', ...
                 'linewidth', 1.5);
            plot(delays(sig), gp_params.tau(sig), ...
                 'marker', '.', ...
                 'linestyle', 'none', ...
                 'color', colors.grays{2}, ...
                 'markersize', pointsize);
            plot(delays_example(sig_example), gp_params_example.tau(sig_example), ...
                 'marker', '^', ...
                 'linestyle', 'none', ...
                 'color', 'k', ...
                 'markerfacecolor', '#D35FBC', ...
                 'markersize', pointsize);

            % Determine which unique latents are shared between the current pair of
            % groups
            sig = find(all(noGlobal([i j],:),1));
            sig_example = find(all(noGlobal_example([i j],:),1));

            subplot(numGroups,numGroups,(i-1)*numGroups+j);
            hold on;
            xlabel(sprintf('Delay from %s to %s (%s)',groupNames{i},groupNames{j},units));
            ylabel(sprintf('GP timescale %s', units));
            xlim([-1.1*maxDelay,1.1*maxDelay]);
            ylim([0,1.1*maxTau]);
            plot(delays(sig), gp_params.tau(sig), ...
                 'marker', '.', ...
                 'linestyle', 'none', ...
                 'color', colors.grays{6}, ...
                 'markersize', pointsize);
             plot(delays_example(sig_example), gp_params_example.tau(sig_example), ...
                 'marker', '^', ...
                 'linestyle', 'none', ...
                 'color', 'k', ...
                 'markerfacecolor', '#D35FBC', ...
                 'markersize', pointsize);
            axis square;
            hold off;
            
        end
    end
end