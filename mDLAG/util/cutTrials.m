function seqOut = cutTrials(seqIn, varargin)
%
% seqOut = cutTrials(seqIn, ...)  
%
% Extracts trial segments that are all of the same length.  Uses
% overlapping segments if trial length is not integer multiple
% of segment length.  Ignores trials with length shorter than 
% one segment length. If no trials are extracted (because the given
% segLength is larger than all raw trial lengths), then this function 
% returns seqIn unaltered.
%
% INPUTS:
%
% seqIn       - data structure, whose nth entry (corresponding to
%               the nth experimental trial) has fields
%                 trialId      -- unique trial identifier
%                 T (1 x 1)    -- number of timesteps in trial
%                 y (yDim x T) -- neural data
%
% OUTPUTS:
%
% seqOut      - data structure, whose nth entry (corresponding to
%               the nth segment) has fields
%                 trialId      -- identifier of trial from which 
%                                 segment was taken
%                 segId        -- segment identifier within trial
%                 T (1 x 1)    -- number of timesteps in segment
%                 y (yDim x T) -- neural data
%
% OPTIONAL ARGUMENTS:
%
% segLength   - length of segments to extract, in number of timesteps.
%               If infinite, entire trials are extracted, i.e., no 
%               segmenting. (default: Inf)
% verbose     - set true to print warning messages. (default: false)
%
% @ 2009 Byron Yu -- byronyu@stanford.edu

  segLength = Inf; % number of timesteps in each segment
  verbose = false;
  assignopts(who, varargin);

  if isinf(segLength)
    seqOut = seqIn;
    return
  end

  seqOut = [];
  for n = 1:length(seqIn)
    T = seqIn(n).T;
    
    % Skip trials that are shorter than segLength
    if T < segLength
      if verbose
        fprintf('Warning: trialId %4d shorter than one segLength...skipping\n',... 
          seqIn(n).trialId); 
      end
      continue
    end
    
    numSeg = ceil(T/segLength);
    
    if numSeg == 1
      cumOL      = 0;
    else
      totalOL    = (segLength*numSeg) - T;
      probs      = ones(1,numSeg-1)/(numSeg-1);
      % mnrnd is very sensitive to sum(probs) being even slightly
      % away from 1 due to floating point round-off.
      probs(end) = 1-sum(probs(1:end-1));
      randOL     = mnrnd(totalOL, probs);
      cumOL      = [0 cumsum(randOL)];
    end

    seg.trialId = seqIn(n).trialId;
    seg.T       = segLength;    
    
    for s = 1:numSeg
      tStart = -cumOL(s) + segLength * (s-1) + 1;
      
      seg.segId   = s;
      seg.y       = seqIn(n).y(:, tStart:(tStart+segLength-1));
            
      seqOut = [seqOut seg];
    end
  end
  
  if isempty(seqOut)
    if verbose
      fprintf('WARNING: no segments extracted for training.  Defaulting to segLength=Inf.\n');
    end
    seqOut = seqIn;
  end
  