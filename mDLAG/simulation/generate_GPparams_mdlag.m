function [tau, eps, D] ...
    = generate_GPparams_mdlag(binWidth, xDim, numGroups, tau_lim, eps_lim, delay_lim)
% [tau, eps, D] ...
%    = generate_GPparams_mdlag(binWidth, xDim, numGroups, tau_lim, eps_lim, delay_lim)
%
% Description: Randomly generate mDLAG GP parameters, within specified
%              constraints.
%
% Arguments:
%
%     Required:
%
%     binWidth  -- float; intended spike count bin width or sample period 
%                  (in units of time). Assume uniform sampling.
%     xDim      -- int; Number of latent states
%     numGroups -- int; Number of groups
%     tau_lim   -- (1 x 2) array; lower- and upper-bounds of GP timescales
%     eps_lim   -- (1 x 2) array; lower- and upper-bounds of GP noise variances
%     delay_lim -- (1 x 2) array; lower- and upper-bounds of delays, in
%                  units of time
%
% Outputs:
%
%     tau       -- (1 x xDim) array; Gaussian process (GP) timescales (ms)
%     eps       -- (1 x xDim) array; GP noise variances
%     D         -- (numGroups x xDim) array; delay matrix (ms)
% 
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     23 Oct 2022 -- Initial full revision.

% Establish limits of each parameter
min_tau = tau_lim(1);
max_tau = tau_lim(2);
min_eps = eps_lim(1);
max_eps = eps_lim(2);
min_delay = delay_lim(1);
max_delay = delay_lim(2);

% GP timescales and noise variances
tau = min_tau + (max_tau-min_tau).*rand(1,xDim); 
% Deal with noise variances on a log scale
eps = exp(log(min_eps) + (log(max_eps)- log(min_eps)).*(rand(1,xDim)));

% Delays
delays = cell(1,numGroups);
for groupIdx = 1:numGroups
    if groupIdx <= 1
        % Set delays to first group to 0 time steps
        delays{groupIdx} = zeros(1,xDim); 
    else
        % All other groups have non-zero delays
        delays{groupIdx} = min_delay + (max_delay-min_delay).*rand(1,xDim);
    end
end
% Fill in output structure
D = cat(1, delays{:});