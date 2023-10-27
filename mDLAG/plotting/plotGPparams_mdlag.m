function gp_params = plotGPparams_mdlag(params,binWidth,varargin)
%
% gp_params = plotGPparams_mdlag(params,binWidth,...)
%
% Description: Plot mDLAG Delay Matrix and latent timescales. Delays and 
%              timescales are plotted as ordered pairs.
%
% Arguments: 
%
%     Required:
%
%     params  -- Structure containing mDLAG model parameters.
%                Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    Cs      -- (1 x numGroups) cell array; List of factor 
%                               loadings 
%                               {(y1Dim x xDim), (y2Dim x xDim), ...}
%                    alphas  -- (numGroups x xDim) array; alpha parameter 
%                               values
%                    phis    -- (1 x numGroups) cell array; List of 
%                               observation precision parameters 
%                               {(y1Dim x 1), (y2Dim x 1), ...}
%                    ds      -- (1 x numGroups) cell array; List of data 
%                               means
%                               {(y1Dim x 1), (y2Dim x 1), ...}
%                    gamma   -- (1 x xDim) array; GP timescales in units of 
%                               time are given by 'binWidth ./ sqrt(gamma)'                                                    
%                    eps     -- (1 x xDim) array; GP noise variances
%                    D       -- (numGroups x xDim) array; delays from 
%                               latents to observed variables. NOTE: Delays
%                               are reported as (real-valued) number of 
%                               time-steps.
%                    xDim    -- int; number of latent variables
%                    yDims   -- (1 x numGroups) array; dimensionalities of 
%                               each observed group
%
%     binWidth   -- float; bin width or sample period (in e.g., ms)
%
%     Optional:
%
%     sigDims  -- (numGroups x xDim) logical array; sigDims(i,j) is 1 if
%                 latent j explains significant shared variance in group i,
%                 0 otherwise. (default: [])
%     units    -- string; units of time of binWidth (for labels)
%                 (default: '')
%
% Outputs:
%     
%    gp_params -- structure containing mDLAG GP parameters, converted into
%                 units of time.
%                 D   -- (numGroups x xDim) array; delays from latents to
%                        observed variables
%                 tau -- (1 x xDim) array; GP timescales
%              
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     04 Oct 2022 -- Initial full revision.

% Set optional arguments
sigDims = [];
units = '';
assignopts(who,varargin);

xDim = params.xDim;
yDims = params.yDims;
numGroups = length(yDims);
if isempty(sigDims)
    % Plot GP parameters for all latents in all groups, regardless of 
    % significance
    sigDims = ones(numGroups,xDim);
end

colors = generateColors(); % Generate custom plotting colors
pointsize = 25; % Size of scatterplot points
% Format units for axis labels
if ~isempty(units)
    units = sprintf(' (%s)', units); 
end

% Convert GP params into units of time
gp_params = getGPparams_mdlag(params, binWidth);
maxTau = max(gp_params.tau);

% Plot shared latents and pairwise delays.
% But don't make a plot if there are no shared latents.
plotShared = false;
for i = 1:numGroups
    for j = 1:numGroups
        % Consider only different groups
        if i ~= j
            % Determine which latents are shared between the current pair
            % of groups
            sig = find(all(sigDims([i j],:),1));
            if ~isempty(sig)
                plotShared = true; 
            end
        end
    end
end

if plotShared
    figure;
    for i = 1:numGroups
        for j = 1:numGroups

            % Consider only different groups
            if i ~= j
                % Take only the delays associated with the current pair of groups
                delays = gp_params.D(j,:) ...
                       - gp_params.D(i,:);

                % Determine which latents are shared between the current pair of
                % groups
                sig = find(all(sigDims([i j],:),1));

                % Figure out plot limits
                maxDelay = max(abs(delays));
                if maxDelay <= 0.05*maxTau
                    maxDelay = 0.5*maxTau;
                end

                subplot(numGroups,numGroups,(i-1)*numGroups+j);
                hold on;
                xlabel(sprintf('Delay from group %d to group %d%s',i,j,units));
                ylabel(sprintf('GP timescale%s', units));
                xlim([-1.1*maxDelay,1.1*maxDelay]);
                ylim([0,1.1*maxTau]);
                line([0 0], [0 1.1*maxTau], ...
                     'Color', colors.grays{6}, ...
                     'linestyle', '--', ...
                     'linewidth', 1.5);
                plot(delays(sig), gp_params.tau(sig), ...
                     'marker', '.', ...
                     'linestyle', 'none', ...
                     'color', colors.grays{1}, ...
                     'markersize', pointsize);
                hold off;
            end
        end
    end
end
% Plot local latents (i.e., latents present in only one group)
% First determine how many local latents are in each group, and
% don't attempt to plot if there are no local latents
localDims = cell(1,numGroups);
plotLocal = false;
for groupIdx = 1:numGroups
    local = zeros(numGroups,1);
    local(groupIdx) = 1;
    localDims{groupIdx} = find(ismember(sigDims',local','rows'));
    if ~isempty(localDims{groupIdx})
        plotLocal = true;
    end
end

if plotLocal
    figure;
    for groupIdx = 1:numGroups
        numLocal = length(localDims{groupIdx});
        % Don't try to plot anything for groups with 0 local latents
        if numLocal > 0
            subplot(1,numGroups,groupIdx);
            hold on;
            xlim([0,numLocal+1]); 
            ylim([0,1.1*maxTau]);
            h = bar(1:numLocal,gp_params.tau(localDims{groupIdx}),0.4);    
            set(h,'facecolor',colors.grays{1},'edgecolor',colors.grays{1});
            ylabel(sprintf('GP timescale%s', units));
            set(gca,'XTick',1:numLocal);
            xlabel(sprintf('Local latents, group %d', groupIdx)); 
            hold off;
        end
    end
end
