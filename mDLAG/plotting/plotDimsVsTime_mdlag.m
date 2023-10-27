function plotDimsVsTime_mdlag(seq, xspec, params, binWidth, varargin)
%
% plotDimsVsTime_mdlag(seq, xspec, binWidth, varargin)
%
% Description: Plot each mDLAG dimension versus time in a separate panel, 
%              along with the mean trajectory across trials.
%              Group and scale panels according to which observation group
%              trajectories belong to.
%
% Arguments:
%
%     Required:
%
%     seq       -- data structure containing extracted trajectories
%     xspec     -- string; field name of trajectories in 'seq' to be 
%                  plotted (e.g., 'xorth' or 'xsm')
%     binWidth  -- bin width used when fitting model
%     params    -- Structure containing mDLAG model parameters.
%                    D       -- (numGroups x xDim) array; delays from 
%                               latents to observed variables. NOTE: Delays
%                               are reported as (real-valued) number of 
%                               time-steps.
%                    xDim    -- int; number of latent variables
%                    yDims   -- (1 x numGroups) array; dimensionalities of 
%                               each observed group
%
%     Optional:
%
%     nPlotMax  -- int; maximum number of trials to plot (default: 20)
%                  NOTE: Not relevant if trialGroups is specified.
%     plotSingle -- logical; if true, plot single-trial trajectories
%                   (default: true)
%     plotMean   -- logical; if true, plot mean across single-trial
%                   trajectories. Mean will be over the nPlotMax trials
%                   being plotted. (default: true)
%     plotErr    -- int; 0: do not plot error bands; 1: plot SEM across 
%                   single-trial trajectories; 2: plot 95% confidence 
%                   bands over single-trial trajectories. Error bands will 
%                   be over the nPlotMax trials being plotted. (default: 0)
%     units     -- string; units of time of binWidth (for labels)
%                  (default: '')
%     trialGroups -- (1 x numTrialGroups) cell array; Each element contains
%                    a list of trials to be grouped together for
%                    color-coding and for computing a mean time course.
%                    (default: {})
%     trialColors -- (1 x numTrialGroups) cell array; If trialGroups is
%                    specified, then trialColors must be specified, 
%                    where each element is the color for a given trial
%                    group. (default: {})
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     04 Oct 2022 -- Initial full revision.

nPlotMax   = 20;
plotSingle = true;
plotMean   = true;
plotErr    = 0;
units      = '';
trialGroups = {};
trialColors = {};
assignopts(who, varargin);

colors = generateColors(); % Get custom colors for plotting
numGroups = length(params.yDims);
xDim = params.xDim;

% Check if there are any trajectories
Xall = [seq.(xspec)];

% Set up figure and axes
if isempty(Xall)
    fprintf('plotDimsVsTime_mdlag: No trajectories to plot.\n');
    return;
end

% Number of trajectories to be plotted
N = 0;
allTrials = [];
if isempty(trialGroups)
    N = min(length(seq), nPlotMax); 
    allTrials = 1:N;
else
    for trialGroupIdx = 1:length(trialGroups)
        N = N + length(trialGroups{trialGroupIdx}); 
        allTrials = [allTrials trialGroups{trialGroupIdx}];
    end
end
  
% Group sequences and parameters
groupSeq = partitionSeq(seq,xDim.*ones(1,numGroups),'datafield',xspec);

% Convert delay labels to units of time
gp_params = getGPparams_mdlag(params, binWidth);
D = gp_params.D;

% Global figure properties
f = figure;

% Number of rows and columns
nRows = numGroups;
nCols = xDim;

% Size and axis scales
set(f, 'Units', 'Normalized', ...
    'OuterPosition', [0.05 0.05 min([1 0.2*nCols]) min([1 0.25*nRows])]);
Tmax    = max([seq.T]);  
Tmin    = min([seq.T]);
xtkStep = ceil(Tmax/25)*5;
xtk     = 1:xtkStep:Tmax;
xtkl    = 0:(xtkStep*binWidth):(Tmax-1)*binWidth;

% Format units for axis labels
if ~isempty(units)
    units = sprintf(' (%s)', units); 
end

fontsize = 12; % Size of text label fonts

for groupIdx = 1:numGroups

    % Determine vertical scale for current group
    Xall = [groupSeq{groupIdx}(allTrials).(xspec)];
    xMax = ceil(10*max(abs(Xall(:)))) / 10;
    xMid = xMax / 2; %ceil(10*(xMax/2)) / 10;
    ytk     = [-xMax -xMid 0 xMid xMax];

    % Across-group latents
    for k = 1:xDim
        h = subplot(nRows, nCols, (groupIdx-1)*nCols+k);
        hold on;
        if isempty(trialGroups)
            dat = cat(3,groupSeq{groupIdx}(allTrials).(xspec));
            for n = 1:length(allTrials)
                if plotSingle 
                    % Plot single-trial trajectories
                    T = groupSeq{groupIdx}(allTrials(n)).T;
                    col = colors.grays{5};
                    plot(1:T, dat(k,:,n), ...
                         'linewidth', 0.05, ...
                         'color', col);
                end
            end
            % Only average over time points up to Tmin, if trial
            % lengths are different.
            xmean = mean(dat(k,1:Tmin,:),3);
            % Plot the mean trajectory
            if plotMean
                plot(1:Tmin, xmean, ...
                     'linewidth', 2.0, ... 
                     'color', colors.grays{1});
            end     
            % Plot error bands
            switch plotErr
                % Only compute error bands over time points up to Tmin, if
                % trial lengths are different.
                case 1
                    % SEM
                    xsem = std(dat(k,1:Tmin,:),0,3) ./ sqrt(length(allTrials));
                    fill([1:Tmin, fliplr(1:Tmin)], [xmean + xsem, fliplr(xmean - xsem)], ...
                         'k', ...
                         'edgecolor', 'none', ...
                         'facealpha', 0.2);
                case 2
                    % 95% confidence bands
                    p = squeeze(prctile(dat(k,1:Tmin,:), [5 95], 3));
                    fill([1:Tmin, fliplr(1:Tmin)], [p(:,2)', fliplr(p(:,1)')], ...
                         'k', ...
                         'edgecolor', 'none', ...
                         'facealpha', 0.2);
            end
        else
            for trialGroupIdx = 1:length(trialGroups)
                dat = cat(3,groupSeq{groupIdx}(trialGroups{trialGroupIdx}).(xspec));
                for n = 1:length(trialGroups{trialGroupIdx})
                    if plotSingle 
                        % Plot single-trial trajectories
                        T = groupSeq{groupIdx}(trialGroups{trialGroupIdx}(n)).T;
                        plot(1:T, dat(k,:,n), ...
                             'linewidth', 0.05, ...
                             'color', trialColors{trialGroupIdx});
                    end
                end
                % Only average over time points up to Tmin, if trial
                % lengths are different.
                xmean = mean(dat(k,1:Tmin,:),3);
                % Plot the mean trajectory
                if plotMean
                    plot(1:Tmin, xmean, ...
                         'linewidth', 2.0, ... 
                         'color', trialColors{trialGroupIdx});
                end
                % Plot error bands
                switch plotErr
                    % Only compute error bands over time points up to Tmin,
                    %  if trial lengths are different.
                    case 1
                        % SEM
                        xsem = std(dat(k,1:Tmin,:),0,3) ./ sqrt(length(trialGroups{trialGroupIdx}));
                        fill([1:Tmin, fliplr(1:Tmin)], [xmean + xsem, fliplr(xmean - xsem)], ...
                             trialColors{trialGroupIdx}, ...
                             'edgecolor', 'none', ...
                             'facealpha', 0.2);
                    case 2
                        % 95% confidence bands
                        p = squeeze(prctile(dat(k,1:Tmin,:), [5 95], 3));
                        fill([1:Tmin, fliplr(1:Tmin)], [p(:,2)', fliplr(p(:,1)')], ...
                             trialColors{trialGroupIdx}, ...
                             'edgecolor', 'none', ...
                             'facealpha', 0.2);
                end
            end
        end
        % Additional formatting of titles and axis labels.
        axis([1 Tmax 1.1*min(ytk) 1.1*max(ytk)]);
        set(h, 'xtick', xtk, 'xticklabel', xtkl);
        set(h, 'ytick', ytk, 'yticklabel', ytk);
        str = sprintf('$${\\mathbf x}_{%d,%d,:}$$',groupIdx,k);
        str = sprintf('%s, $$D_{%d,%d} = %3.1f$$', str, groupIdx, k, D(groupIdx,k));
        str = sprintf('%s%s', str, units);
        title(str, 'interpreter', 'latex', 'fontsize', fontsize);
        xlabel(sprintf('Time%s', units));
    end

end