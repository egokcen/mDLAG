function params = generate_params_mdlag(yDims, xDim, binWidth, ...
                                        hyperparams, snr, tau, eps, D)
% params = generate_params_mdlag(...)
%
% Description: Randomly generate mDLAG model parameters, within 
%              specified constraints. The params output structure is
%              compatible with other mDLAG codepack functions.
%
% Arguments:
%
%     yDims     -- (1 x numGroups) array; List of dimensionalities of
%                  observed data, [y1Dim, y2Dim, ...]
%     xDim      -- int; Dimensionality of latents, X
%     binWidth  -- float; intended spike count bin width or sample period 
%                  (in units of time). Assume uniform sampling.
%     hyperparams -- structure with the following fields:
%                      beta  -- positive float; precision of mean parameter
%                               generative model (Gaussian)
%                      a_phi -- positive float; 'a' shape parameter of
%                               observation precision (phi) generative 
%                               model (Gamma with mean a/b)
%                      b_phi -- positive float; 'b' scale parameter of
%                               observation precision (phi) generative
%                               model (Gamma with mean a/b)
%                      a_alpha -- (numGroups x xDim) array; 'a' shape 
%                                 parameter of alpha parameter generative
%                                 model (Gamma with mean a/b)
%                      b_alpha -- (numGroups x xDim) array; 'b' scale 
%                                 parameter of alpha parameter generative 
%                                 model (Gamma with mean a/b)
%     snr       -- (1 x numGroups) array; List of signal-to-noise ratios,
%                  defined as trace(CC') / sum(1./phi)
%     tau       -- (1 x xDim) array; Gaussian process (GP) timescales (ms)
%     eps       -- (1 x xDim) array; GP noise variances
%     D         -- (numGroups x xDim) array; delay matrix (ms)
%
% Outputs:
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
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     27 Sep 2022 -- Initial full revision.

numGroups = length(yDims);

% Initialize output structure
params = struct('covType', 'rbf', ...
                'Cs', [], ...
                'alphas', nan(numGroups,xDim), ...
                'phis', [], ...
                'ds', [], ...
                'gamma', (binWidth ./ tau).^2, ...
                'eps', eps, ...
                'D', D./binWidth, ...
                'xDim', xDim, ...
                'yDims', yDims ...
               );
params.Cs = cell(1,numGroups);
params.phis = cell(1,numGroups);
params.ds = cell(1,numGroups);

% Generate observation model parameters
for groupIdx = 1:numGroups
    
    for xIdx = 1:xDim
        params.alphas(groupIdx,xIdx) = gamrnd(hyperparams.a_alpha(groupIdx,xIdx), 1./hyperparams.b_alpha(groupIdx,xIdx));
    end
    params.ds{groupIdx} = mvnrnd(zeros(1,yDims(groupIdx)), hyperparams.beta^(-1).*eye(yDims(groupIdx)), 1)';
    params.phis{groupIdx} = gamrnd(hyperparams.a_phi, 1./hyperparams.b_phi, yDims(groupIdx), 1);
    params.Cs{groupIdx} = nan(yDims(groupIdx),xDim);
    for xIdx = 1:xDim
        params.Cs{groupIdx}(:,xIdx) = mvnrnd(zeros(1,yDims(groupIdx)), params.alphas(groupIdx,xIdx)^(-1).*eye(yDims(groupIdx)), 1)';
    end
    
    % Enforce the desired signal-to-noise ratios
    varCC = trace(params.Cs{groupIdx} * params.Cs{groupIdx}');
    varNoise_desired = varCC / snr(groupIdx);
    varNoise_current = sum(params.phis{groupIdx}.^(-1));
    params.phis{groupIdx} = params.phis{groupIdx} .* (varNoise_current / varNoise_desired);
 
end