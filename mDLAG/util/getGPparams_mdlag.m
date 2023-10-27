function gp_params = getGPparams_mdlag(params, binWidth)
%
% gp_params = getGPparams_mdlag(params, binWidth)
%
% Description: Delays and GP timescales can be found in params, but they
%              are easier to interpret when given in units of time. This
%              function gets GP timescales and delays from params and 
%              converts them into the units of time corresponding to 
%              binWidth. binWidth should match the binWidth to which the 
%              model was fitted.
%
% Arguments: 
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

% Convert delays from bins to units of time
gp_params.D = binWidth .* params.D;

% Convert timescales to units of time
gp_params.tau = binWidth ./ sqrt(params.gamma);
