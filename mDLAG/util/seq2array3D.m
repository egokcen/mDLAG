function data = seq2array3D(seq, varargin)
%
% data = seq2array3D(seq, ...)
%
% Description: Convert seq (mDLAG-compatible) structure to a 3D array.
%
% Arguments:
%
%     Required: 
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
%
%     NOTE: T is assumed to be the same for each trial in seq.
%
%     Optional:
%     
%     datafield -- string; Name of data field in seq (default: 'y')
%
% Outputs:
%
%     data     -- (yDim x T x N) array with activity of each unit 
%                 across time and trials


datafield = 'y';
extraOpts = assignopts(who, varargin);

N = length(seq);
[yDim, T] = size(seq(1).(datafield));

data = nan(yDim, T, N);
for n = 1:N
    data(:,:,n) = seq(n).(datafield);
end