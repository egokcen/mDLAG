function [pairDims, pairVarExp, pairs] = computeDims_pairs(dims, dimTypes, varExp)
%
% [pairDims, pairVarExp, pairs] = computeDims_pairs(dims, dimTypes, varExp)
%
% Description: Compute the total dimensionality of each group and the 
%              shared dimensionality between each pair of groups given
%              the dimensionalities and types of dimensions in dims and
%              dimTypes, respectively. Compute also the shared variance 
%              explained by a pairwise interaction in each area. 
% 
% Arguments:
%
%     dims     -- (1 x numDimTypes) array; The number of each type of
%                 dimension. dims(i) corresponds to the dimension type in 
%                 dimTypes(:,i).
%     dimTypes -- (numGroups x numDimTypes) array; dimTypes(:,i) is a
%                 binary vector indicating the structure of dimension type
%                 i. 1 indicates a group is involved, 0 otherwise.
%     varExp   -- (numGroups x numDimTypes) array; varExp(i,j) is the
%                 proportion of shared variance explained in group i by
%                 dimension type j. varExp(:,j) corresponds to the
%                 dimension type in dimTypes(:,j).
%
% Outputs:
%
%     pairDims   -- (numPairs x 3) array; 
%                       pairDims(i,1): total dimensionality of group 1 in 
%                                      pair i.
%                       pairDims(i,2): shared dimensionality between pair i
%                       pairDims(i,3): total dimensionality of group 2 in 
%                                      pair i.
%     pairVarExp -- (numPairs x 2) array;
%                       pairVarExp(i,1): shared variance explained by
%                                        pairwise interaction, group 1 of
%                                        pair i
%                       pairVarExp(i,2): shared variance explained by
%                                        pairwise interaction, group 2 of
%                                        pair i
%     pairs      -- (numPairs x 2) array; pairs(i,:) gives the indexes of
%                   groups in pair i, and corresponds to the relevant rows
%                   of pairDims and pairVarExp.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     07 Sep 2022 -- Initial full revision.

numGroups = size(dimTypes,1);

% For each group, create a list of all dimension types that involve that
% group
groupIdxs = cell(1,numGroups);
for groupIdx = 1:numGroups
    groupIdxs{groupIdx} = find(dimTypes(groupIdx,:));
end

% Create a list of all possible pairs
pairs = nchoosek(1:numGroups,2);
numPairs = size(pairs,1);

% Compute means
pairDims = nan(numPairs,3);
pairVarExp = nan(numPairs,2);
for pairIdx = 1:numPairs
    % Total in group 1
    pairDims(pairIdx,1) = sum(dims(groupIdxs{pairs(pairIdx,1)}));

    % Shared between both groups
    sharedIdxs = intersect(groupIdxs{pairs(pairIdx,1)}, groupIdxs{pairs(pairIdx,2)});
    pairDims(pairIdx,2) = sum(dims(sharedIdxs));

    % Total in group 2
    pairDims(pairIdx,3) = sum(dims(groupIdxs{pairs(pairIdx,2)}));

    % Variance explained in each group by pairwise interaction
    pairVarExp(pairIdx,1) = sum(varExp(pairs(pairIdx,1),sharedIdxs),2);
    pairVarExp(pairIdx,2) = sum(varExp(pairs(pairIdx,2),sharedIdxs),2);

end