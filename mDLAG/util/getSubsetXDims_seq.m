function seqSub = getSubsetXDims_seq(seq, xDim, numGroups, xDims_kept, varargin)
%
% seqSub = getSubsetXDims_seq(seq, xDim, numGroups, xDims_kept, varargin)
%
% Description: Get sequences corresponding to the dimensions in 
%              'xDims_kept'.
%
% Arguments:
%
%     Required:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                   trialId                -- unique trial identifier
%                   T (1 x 1)              -- number of timesteps
%                   (datafield) (xDim x T) -- latent time courses
%     xDim        -- int; original latent dimensionality, prior to removing
%                    the desired selection
%     numGroups   -- int; number of observation groups
%     xDims_kept  -- (1 x numDims) array; latent state dimensions to be  
%                    retained in seqSub.
%
%     Optional:
%     
%     datafield -- string; Name of data field in seq (default: 'xsm')
%
% Outputs:
%
%     seqSub -- Same structure as 'seq' above, but with only sequences
%               of the specified 'xDims' remaining.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     21 Jan 2022 -- Initial full revision.

datafield = 'xsm';
extraOpts = assignopts(who, varargin);

% Initialize other relevant variables
block_idxs = get_block_idxs(repmat(xDim,1,numGroups));

% Get the sequences corresponding to the latents in 'xDims_kept'
for n = 1:length(seq)
    seqSub(n).trialId = seq(n).trialId;
    seqSub(n).T = seq(n).T;
    seqSub(n).(datafield) = [];
    for groupIdx = 1:numGroups
        currBlockIdxs = block_idxs{groupIdx}(1):block_idxs{groupIdx}(2);
        seqSub(n).(datafield) = [seqSub(n).(datafield); seq(n).(datafield)(currBlockIdxs(xDims_kept),:)]; 
    end
end