function [Range] = compressToRange(IntegerArray)
%compressToRange compresses the 'IntegerArray' to a range of integers.
% This function compresses an array of integers into a range of integers
% with no missing values.

% EXAMPLES:
%   Range = smi_helpers.compressToRange([2; 4; 6; 6; 7; 10])
%       -> Range = [1; 2; 3; 3; 4; 5]
%   Range = smi_helpers.compressToRange([1; 5; 8; 8; 11; 3])
%       -> Range = [1; 3; 4; 4; 5; 2]

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Ensure Range consists of the set of integers 1:numel(IntegerArray)
% without skipping any integers.
[IntsSorted, SortIndices] = sort(IntegerArray);
IntList = unique(IntsSorted);
NPerInt = groupcounts(IntsSorted);
NCumulative = [0; cumsum(NPerInt)];
NInts = numel(IntList);
for nn = 1:NInts
    % Set all entries with IntegerArray==nn equal to nn.
    IndexArray = (1:NPerInt(nn)) + NCumulative(nn);
    IntsSorted(IndexArray) = nn;
end
Range = NaN(size(IntegerArray));
Range(SortIndices) = IntsSorted;


end