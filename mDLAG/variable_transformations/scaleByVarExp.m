function [seq, sortParams] = scaleByVarExp(seq, params, varexp, varargin)
%
% [seq, sortParams] = scaleByVarExp(seq, params, varexp, varargin)
%
% Description: Scale (and re-order) latent trajectories by their fraction 
%              of variance explained, for visualization purposes. 
%
% Arguments:
%
%     Required:
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId        -- unique trial identifier
%                     T (1 x 1)      -- number of timesteps
%                     (indatafield) (xDim x T) -- latent trajectories
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
%     varexp   -- (numGroups x xDim) array; varexp(i,:) is an
%                 array with the variance explained in group i by each
%                 individual latent variable.
%
%     Optional:
%
%     sortDims     -- logical; set to true to sort dimensions by variance
%                     explained. (default: true)
%     sortGroup    -- int; if sorting, choose which group to use as an
%                     anchor when sorting across-group dimensions.
%                     (default: 1)
%     numDim       -- int; Keep at most the top numDim dimensions 
%                     (default: xDim)
%     indatafield  -- string; fieldname of input data (default: 'xsm')
%     outdatafield -- string; fieldname of output data (default: 'xve')
%
% Outputs:
%
%     seq     -- input data structure with new field
%                    (outdatafield) -- (xDim x T) array; scaled (and
%                                      re-ordered) latent trajectories
%     sortParams -- same structure as params, with dimensions re-sorted.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     04 Oct 2022 -- Initial full revision.

sortDims = true;
sortGroup = 1;
xDim = params.xDim;
numDim = xDim;
indatafield = 'xsm';
outdatafield = 'xve';
assignopts(who, varargin);

numGroups = length(params.yDims);
N = length(seq);

% Initialize output structures
sortParams = params;

% Scale latents by variance explained
varexp_all = reshape(varexp',[],1);
for n = 1:N
    seq(n).(outdatafield) = zeros(size(seq(n).(indatafield)));
    seq(n).(outdatafield) = varexp_all .* seq(n).(indatafield);
end

if sortDims
    % Initialize sort indices
    sortIdxs = 1:xDim;
    
    % Sort dimensions by variance explained.
    
    % Sort all dimensions according to the anchor group 
    [~, sortIdxs] = sort(varexp(sortGroup,1:xDim),'descend');
    % Keep at most the top numDim dimensions
    keptDim = min([numDim xDim]);
    sortIdxs = sortIdxs(1:keptDim);
    
    % Sort parameters
    sortParams = getSubsetXDims_params(params, sortIdxs);
    
    % Sort time courses
    groupSeq = partitionSeq(seq,xDim.*ones(1,numGroups),'datafield',outdatafield);
    for n = 1:N
        currSeq = [];
        for groupIdx = 1:numGroups
            currSeq = [currSeq; groupSeq{groupIdx}(n).(outdatafield)(sortIdxs,:)];
        end
        seq(n).(outdatafield) = currSeq;
    end
end