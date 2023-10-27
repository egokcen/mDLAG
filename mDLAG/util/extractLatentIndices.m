function idxs = extractLatentIndices(xDim,T,numGroups)
%
% idxs = extractLatentIndices(xDim,T,numGroups)
%
% Description: Extract indices into K_big or SigX that correspond to
%              individual latents.
%
% Arguments:
%
%     xDim      -- int; number of latent states
%     T         -- int; number of timesteps
%     numGroups -- int; number of groups
%
% Outputs:
%
%     idxs      -- (1 x xDim) cell array; 
%                    idxs(i) -- (1 x numGroups*T) array; indices 
%                               corresponding to latent state i
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     19 Oct 2022 -- Initial full revision.

idxs = cell(1,xDim);
for i = 1:xDim
    idxs{i} = [];
    for t = 1:T
        bigBaseIdx = (xDim*numGroups)*(t-1) + 1;                
        for groupIdx = 1:numGroups
            bigIdx = bigBaseIdx + (groupIdx-1) * xDim;
            idxs{i} = [idxs{i} bigIdx+i-1];
        end
    end
end