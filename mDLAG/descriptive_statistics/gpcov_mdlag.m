function k = gpcov_mdlag(params, binWidth, maxlag, varargin)
%
% k = gpcov_mdlag(params, binWidth, maxlag, ...)
%
% Description: Compute and (optionally) plot GP covariance functions over
%              a range of time lags specified by maxlag.
%
% Arguments:
%
%     Required:
%
%     params  -- Structure containing mDLAG model parameters.
%                Contains the fields
%         covType    -- string; type of GP covariance (e.g., 'rbf')
%         gamma      -- (1 x xDim) array; GP timescales in units of 
%                       time are given by 'binWidth ./ sqrt(gamma)'                                                    
%         eps        -- (1 x xDim) array; GP noise variances
%         D          -- (numGroups x xDim) array; delays from latents to 
%                       observed variables. NOTE: Delays are reported as 
%                       (real-valued) number of time-steps.
%         d.mean     -- (yDim x 1) array; posterior mean of mean parameter
%         d.cov      -- (yDim x 1) array; diagonal elements of the
%                       posterior covariance matrix of d
%         C.means    -- (numGroups x 1) cell array; yDims(groupIdx) x xDim
%                       mean loadings matrix for each group
%         C.covs     -- (numGroups x 1) cell array; C.covs{groupIdx) is a
%                       (yDims(groupIdx) x 1) cell array, and each element
%                       is a (xDim x xDim) matrix giving the posterior
%                       covariance of a row of C.
%         C.moments  -- (numGroups x 1) cell array; C.moments{groupIdx) is a
%                       (yDims(groupIdx) x 1) cell array, and each element
%                       is a (xDim x xDim) matrix giving the posterior
%                       second moment of a row of C.
%         alpha.a    -- (numGroups x 1) array; shape parameters of alpha 
%                       posterior
%         alpha.b    -- (numGroups x xDim) array; scale parameters of 
%                       alpha posterior
%         alpha.mean -- (numGroups x xDim) array; mean precisions of
%                       loading weights (for ARD); alpha.a ./ alpha.b
%         phi.a      -- float; shape parameter of phi posterior
%         phi.b      -- (yDim x 1) array; scale parameters of phi posterior
%         phi.mean   -- (yDim x 1) array; mean precisions of observations; 
%                       alpha.a ./ alpha.b
%         xDim       -- int; number of latent variables
%         yDims      -- (1 x numGroups) array; dimensionalities of each 
%                       observed group
%
%     binWidth   -- float; resolution (sample period or bin width), in
%                   units of time, at which params were estimated.
%     maxlag     -- float; maximum time lag to consider when computing
%                   covariance function; same units as binWidth.
%
%     Optional:
%
%     showPlot -- logical; set true to plot all transformed covariance
%                 functions (default: true)
%     stepres  -- float; resolution of the computed covariance 
%                 functions, in the same units of time as binWidth
%                 (default: 0.1).
%
% Outputs:
%
%     k -- (numGroups x numGroups) cell array; k{i,j} is a 
%          (1 x xDim_across) cell array, and k{i,j}{r} is the cross-
%          covariance function between group i and group j, given by
%          latent r.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     01 Feb 2023 -- Initial full revision.

showPlot = true;
stepres = 0.1;
assignopts(who,varargin);

numGroups = length(params.yDims);
xDim = params.xDim;
lagSteps = -maxlag:stepres:maxlag;
T = length(lagSteps);

% Convert GP params to same units as binWidth
gp_params = getGPparams_mdlag(params,binWidth);

% Compute prior covariance functions
k = cell(numGroups,numGroups);
for groupIdx1 = 1:numGroups
    for groupIdx2 = 1:numGroups
        k{groupIdx1,groupIdx2} = cell(1,xDim);
        for xIdx = 1:xDim
            D = (gp_params.D(groupIdx2,xIdx) - gp_params.D(groupIdx1,xIdx));
            k{groupIdx1,groupIdx2}{xIdx} ...
                = (1 - params.eps(xIdx)).*exp(-0.5*(gp_params.tau(xIdx)^(-2)).*(lagSteps-D).^2);
            k{groupIdx1,groupIdx2}{xIdx}((lagSteps - D) == 0) ...
                = k{groupIdx1,groupIdx2}{xIdx}((lagSteps - D) == 0) ...
                  + params.eps(xIdx);
        end
    end
end

if showPlot
    % Across-group covariance functions
    for xIdx = 1:xDim
        figure;
        for groupIdx1 = 1:numGroups
            for groupIdx2 = 1:numGroups
                k_curr = k{groupIdx1,groupIdx2}{xIdx};
                subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                hold on;
                plot(lagSteps, k_curr, 'k-', 'linewidth', 1.5);
                [maxval, maxdelay] = max(k_curr);
                line([lagSteps(maxdelay) lagSteps(maxdelay)], [0 1.05*max([1 maxval])], 'color', 'r', 'linestyle', '--');
                line([0 0], [0 1.05*max([1 maxval])], 'color', 'k', 'linestyle', '--');
                axis square;
                xlabel('Time lag (ms)');
                ylabel('Correlation');
                axis([min(lagSteps) max(lagSteps) 0 1.05*max([1 maxval])]);
            end
        end
    end
end