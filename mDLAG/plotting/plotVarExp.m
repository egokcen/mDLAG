function plotVarExp(varExp, dimTypes, varargin)
%
% plotVarExp(varExp, dimTypes, ...)
%
% Description: Plot the shared variance explained by each type of dimension
%             (given by dimTypes) in each group.
% 
% Arguments:
%
%     Required:
%
%     varExp   -- (numGroups x numDimTypes) array; varExp(i,j) is the
%                 proportion of shared variance explained in group i by
%                 dimension type j. varExp(:,j) corresponds to the
%                 dimension type in dimTypes(:,j).
%     dimTypes -- (numGroups x numDimTypes) array; dimTypes(:,i) is a
%                 binary vector indicating the structure of dimension type
%                 i. 1 indicates a group is involved, 0 otherwise.
%
%     Optional:
%
%     varExp_sem  -- (numGroups x numDimTypes) array; The standard error of  
%                    each element in varExp, if varExp is the mean over 
%                    different runs. (default: [])
%     groupNames  -- (1 x numGroups) cell array; List of strings containing 
%                    the name of each group. By default, names will be '1', 
%                    '2', '3', etc.
%     plotZeroDim -- logical; Set true to plot the number of dimensions
%                    that are not significant in any group. (default:
%                    false)
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

varExp_sem = [];
[numGroups, numDimTypes] = size(dimTypes);
groupNames = cell(1,numGroups);
for groupIdx = 1:numGroups
    groupNames{groupIdx} = num2str(groupIdx); 
end
plotZeroDim = false;
assignopts(who, varargin);

% Sort dimension types by increasing number of involved groups
dimCard = sum(dimTypes,1); % Cardinality of each dimension type
dimSortIdxs = [];
if plotZeroDim
    % Visualize dimensions significant in no groups
    groupList = 0:numGroups;
else
    % Skip dimensions significant in no groups
    groupList = 1:numGroups;
end

for groupIdx = groupList
    currDims = find(dimCard == groupIdx);
    dimSortIdxs = [dimSortIdxs currDims];
end

% Set up labels for x-axis
xtcklbls = cell(1,numDimTypes);
for dimIdx = 1:numDimTypes
    if dimIdx <= 1
        xtcklbls{dimIdx} = 'n.s.';  % 'Not significant'
    else
        xtcklbls{dimIdx} = '';
        for groupIdx = 1:numGroups
            if dimTypes(groupIdx,dimIdx)
               xtcklbls{dimIdx} = [xtcklbls{dimIdx} groupNames{groupIdx} '-'];
            end
        end
        xtcklbls{dimIdx} = xtcklbls{dimIdx}(1:end-1);
    end
end

% Plot variance explained in each group
figure;
for groupIdx = 1:numGroups
    subplot(numGroups,1,groupIdx);
    hold on;
    bar(varExp(groupIdx,dimSortIdxs), ...
        'facealpha', 0.3, ...
        'facecolor', 'k', ...
        'edgecolor', 'k', ...
        'linewidth', 1.5);
    if ~isempty(varExp_sem)
        errorbar(1:length(dimSortIdxs), varExp(groupIdx,dimSortIdxs), varExp_sem(groupIdx,dimSortIdxs), ...
             'k.', ...
             'linewidth', 1.5, ...
             'marker', 'none');
    end
    xlabel('Dimension type');
    ylabel('Frac. Shared Var. Exp.');
    xticks(1:length(dimSortIdxs));
    xticklabels(xtcklbls(dimSortIdxs));
    title(groupNames{groupIdx});
    ylim([0 1]);
end
