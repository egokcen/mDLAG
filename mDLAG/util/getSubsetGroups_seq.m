function seqSub = getSubsetGroups_seq(seq, dims, groups, varargin)
%
% seqSub = getSubsetGroups_seq(seq, dims, groups, ...)
%
% Description: Get sequences corresponding to the subset of groups in 
%              'groups'.
%
% Arguments:
%
%     Required:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                   trialId               -- unique trial identifier
%                   T (1 x 1)             -- number of timesteps
%                   (datafield) (dim x T) -- continuous valued data
%     dims   -- (1 x numGroups) array; the dimensionalities of each 
%               original group in seq
%     groups -- (1 x numKeptGroups) array; list of indices of the desired
%               group subset.
%
%     Optional:
%     
%     datafield -- string; Name of data field in seq (default: 'y')
%
% Outputs:
%
%     seqSub -- Same structure as 'seq' above, but with only sequences
%               in the specified 'groups' remaining.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     22 Oct 2022 -- Initial full revision.

datafield = 'y';
extraOpts = assignopts(who, varargin);

% Initialize other relevant variables
block_idxs = get_block_idxs(dims);

% Get all indices belonging to the desired groups
keptGroupIdxs = [];
for groupIdx = groups
    currGroup = block_idxs{groupIdx};
    keptGroupIdxs = [keptGroupIdxs currGroup(1):currGroup(2)];
end

% Get the sequences corresponding to the groups in 'groups'
for n = 1:length(seq)
    seqSub(n).trialId = seq(n).trialId;
    seqSub(n).T = seq(n).T;
    seqSub(n).(datafield) = seq(n).(datafield)(keptGroupIdxs,:); 
end