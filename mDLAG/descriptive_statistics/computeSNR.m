function snr = computeSNR(C, phi)
%
% snr = computeSNR(C, phi)
%
% Description: Compute the signal-to-noise ratio (SNR) of each observed
%              group, according to mDLAG model parameters.
% 
% Arguments:
%
%     C.moments  -- (numGroups x 1) cell array; C.moments{groupIdx) is a
%                   (yDims(groupIdx) x 1) cell array, and each element
%                   is a (xDim x xDim) matrix giving the posterior
%                   second moment of a row of C.
%     phi.mean   -- (yDim x 1) array; mean precisions of observations
%
% Outputs:
%
%     snr -- (numGroups x 1) array; The SNR of each observed group.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     27 Sep 2022 -- Initial full revision.

numGroups = length(C.moments);
yDims = nan(1,numGroups);
for groupIdx = 1:numGroups
    yDims(groupIdx) = length(C.moments{groupIdx}); 
end
obsIdxs = get_block_idxs(yDims);

snr = nan(numGroups,1);
for groupIdx = 1:numGroups
    CC_m = zeros(size(C.moments{groupIdx}{1}));
    for yIdx = 1:yDims(groupIdx)
        CC_m = CC_m + C.moments{groupIdx}{yIdx}; 
    end
    obsBlock = obsIdxs{groupIdx};
    snr(groupIdx) = trace(CC_m) / sum(1./phi.mean(obsBlock(1):obsBlock(2)));
end