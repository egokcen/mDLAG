function plotDims_pairs(pairDims, pairs, numGroups, varargin)
%
% plotDims_pairs(pairDims, pairs, numGroups, ...)
%
% Description: Visualize pairwise analyses of dimensionality: the total
%              dimensionality of each group in each pair, and the shared 
%              dimensionality between each pair.
% 
% Arguments:
%
%     Required:
%
%     pairDims   -- (numPairs x 3) array; 
%                       pairDims(i,1): total dimensionality of group 1 in 
%                                      pair i.
%                       pairDims(i,2): shared dimensionality between pair i
%                       pairDims(i,3): total dimensionality of group 2 in 
%                                      pair i.
%     pairs      -- (numPairs x 2) array; pairs(i,:) gives the indexes of
%                   groups in pair i, and corresponds to the relevant rows
%                   of pairDims and pairVarExp.
%     numGroups  -- int; number of observed groups
%
%     Optional:
%
%     pairDims_sem -- (numPairs x 3) array; The standard error of each 
%                     element in pairDims, if pairDims is the mean over 
%                     different runs. (default: [])
%     groupNames   -- (1 x numGroups) cell array; List of strings containing 
%                     the name of each group. By default, names will be '1', 
%                     '2', '3', etc.
%
% Outputs:
%
%     none.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     07 Sep 2022 -- Initial full revision.
%     22 Sep 2022 -- Added option to plot SEM.

pairDims_sem = [];
groupNames = cell(1,numGroups);
for groupIdx = 1:numGroups
    groupNames{groupIdx} = num2str(groupIdx); 
end
assignopts(who, varargin);

numPairs = size(pairs,1);

% Set up x-axis labels
xtcklbls = cell(numPairs,3);
for pairIdx = 1:numPairs
    % Total in group 1
    xtcklbls{pairIdx,1} = sprintf('Total, %s', groupNames{pairs(pairIdx,1)});
    % Shared between both groups
    xtcklbls{pairIdx,2} = sprintf('%s-%s', groupNames{pairs(pairIdx,1)}, groupNames{pairs(pairIdx,2)});
    % Total in group 2
    xtcklbls{pairIdx,3} = sprintf('Total, %s', groupNames{pairs(pairIdx,2)});
end

% Plot dimensionalities
figure;
for pairIdx = 1:numPairs
    subplot(1,numPairs,pairIdx);
    hold on;
    bar(pairDims(pairIdx,:), ...
        'facealpha', 0.3, ...
        'facecolor', 'k', ...
        'edgecolor', 'k', ...
        'linewidth', 1.5);
    if ~isempty(pairDims_sem)
        errorbar(1:length(pairDims(pairIdx,:)), pairDims(pairIdx,:), pairDims_sem(pairIdx,:), ...
                 'k.', ...
                 'linewidth', 1.5, ...
                 'marker', 'none');
    end
    xlabel('Dimension type');
    ylabel('Dimensionality');
    xticks(1:length(pairDims(pairIdx,:)));
    xticklabels(xtcklbls(pairIdx,:));
    title(xtcklbls{pairIdx,2});
end
