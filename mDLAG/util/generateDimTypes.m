function dimTypes = generateDimTypes(numGroups)
%
% dimTypes = generateDimTypes(numGroups)
%
% Description: Helper function to generate an array with all types of
%              dimensions (singlet, doublet, triplet, global, etc.) given
%              the number of groups, numGroups.
% 
% Arguments:
%
%     numGroups -- int; number of observed groups
%
% Outputs:
%
%     dimTypes -- (numGroups x numDimTypes) array; dimTypes(:,i) is a
%                 binary vector indicating the structure of dimension type
%                 i. 1 indicates a group is involved, 0 otherwise.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     07 Sep 2022 -- Initial full revision.

numDimTypes = 2^numGroups; % Number of dimension types
dimTypes = nan(numGroups,numDimTypes);
for dimIdx = 0:numDimTypes-1
    dimStr = dec2bin(dimIdx,numGroups);
    dimType = nan(numGroups,1);
    for groupIdx = 1:numGroups
        dimType(groupIdx) = str2num(dimStr(groupIdx));
    end
    dimTypes(:,dimIdx+1) = dimType;
end