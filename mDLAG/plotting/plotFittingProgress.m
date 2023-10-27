function plotFittingProgress(trackedParams, binWidth, varargin)
%
% plotFittingProgress(trackedParams, binWidth, varargin)
%
% Description: Plot intermediate values stored throughout mDLAG model
%              fitting vs iteration number. Visualizing these values can
%              help with troubleshooting.
%
% Arguments:
%
%     Required:
%
%     trackedParams -- structure containing parameters tracked  
%                      throughout fitting:
%         lb       -- (1 x numIters) array; variational lower bound 
%                      at each iteration.
%                      Note: Entries will be NaN, on iterations
%                      was lb was not computed, to save time.
%         iterTime -- (1 x numIters) array; computation time for 
%                     each EM iteration
%         Ds       -- (1 x numIters) cell array; the estimated delay 
%                     matrix (D) after each EM iteration.
%         gams     -- (1 x numIters) cell arry; estimated gamma  
%                     after each EM iteration.
%         alphas   -- (1 x numIters) cell arry; estimated ARD 
%                     parameters (alpha.mean) after each EM
%                     iteration.
%     binWidth -- float; bin width or sample period, in units of time 
%                 (e.g., ms)
%
%     Optional:
%     
%     freqLB    -- int; if known, specify how often LL was computed during
%                  model fitting (default: 1)
%     freqParam -- int; if known, specify how often GP parameters were
%                  stored during model fitting (default: 1)
%     units     -- string; units of time of binWidth (for labels)
%                  (default: '')
%
% Outputs:
%     None. (But creates a figure)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     21 Oct 2022 -- Initial full revision.
%     08 Dec 2022 -- Changed input argument structure.

freqLB = 1;
freqParam = 1;
units = '';
assignopts(who,varargin);

% Format units for axis labels
if ~isempty(units)
    units = sprintf(' (%s)', units); 
end

[numGroups,xDim] = size(trackedParams.Ds{1});

colors = generateColors(); % Get custom colors for plotting
fontsize = 12;

figure;

% Lower bound curve
subplot(1,2,1);
hold on;
lb = trackedParams.lb(~isnan(trackedParams.lb));
xlb = [];
if freqLB <= 1
    xlb = 1:length(lb);
elseif freqLB == 2
    xlb = [1 freqLB.*(1:length(lb)-1)];
else
    xlb = [1 2 freqLB.*(1:length(lb)-2)];
end
plot(xlb, lb, 'color', colors.grays{1}, 'linewidth', 1.5);
xlabel('# Iterations');
ylabel('Lower bound');

% Cumulative fitting time
subplot(1,2,2);
hold on;
plot(cumsum(trackedParams.iterTime), 'color', colors.grays{1}, 'linewidth', 1.5);
xlabel('# Iterations');
ylabel('Cumulative time elapsed (s)');

if xDim > 0
    % GP timescale progress
    figure; 
    hold on;
    Gams = cat(3, trackedParams.gams{:});
    for i = 1:xDim
        plot(freqParam.*(0:length(trackedParams.gams)-1), ...
             binWidth./sqrt(squeeze(Gams(1,i,:))), ...
             'linewidth', 1.5);
    end
    xlabel('# Iterations');
    ylabel(sprintf('GP timescale%s', units));
    
    if numGroups > 1
        % Delay progress
        figure;
        Delays = cat(3, trackedParams.Ds{:});
        for rowIdx = 1:numGroups
            for colIdx = 1:numGroups
                % Avoid redundant plots, and plot only different groups
                if colIdx > rowIdx
                    subplot(numGroups,numGroups,numGroups*(rowIdx-1)+colIdx);
                    hold on;
                    for i = 1:xDim
                        plot(freqParam.*(0:length(trackedParams.Ds)-1), ...
                             binWidth.*squeeze(Delays(colIdx,i,:) - Delays(rowIdx,i,:)),...
                             'linewidth', 1.5);
                    end
                    xlabel('# Iterations');
                    ylabel(sprintf('Delay%s, group %d to %d', units, rowIdx, colIdx));
                end
            end
        end
    end
end

% ARD parameter progress
alphas = cat(3, trackedParams.alphas{:});
figure;
for groupIdx = 1:numGroups
    subplot(1,numGroups,groupIdx);
    hold on;
    for i = 1:xDim
        plot(freqParam.*(0:length(trackedParams.alphas)-1), ...
             1./squeeze(alphas(groupIdx,i,:)), ...
             'linewidth', 1.5);
    end
    xlabel('# Iterations');
    ylabel('\alpha^{-1}', 'fontsize', fontsize);
    title(sprintf('Group %d', groupIdx));
end