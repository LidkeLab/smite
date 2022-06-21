function [StartInds, EndInds] = findStartEndInds(BoolArray)
%findStartEndInds finds the start and end indices of events in BoolArray.
% This function finds the starting and ending indices of each series of
% trues in a 1D boolean array.  For example,
% findStartEndInds([0; 1; 1; 1; 0; 0; 1; 0; 1; 1]) should return
% StartInds = [2; 7; 9] and EndInds = [4; 7; 10].
%
% INPUTS:
%   BoolArray: 1D array of logical values.
%
% OUTPUTS:
%   StartInds: Array of starting indices for series of trues.
%   EndInds: Array of ending indices corresponding to the start indices in
%            'StartInds'.

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Find the starting indices of each series of trues, ensuring that we
% capture events starting with the first index of BoolArray.
StartInds = unique([...
    find(BoolArray, 1, 'first'); ...
    find(diff(BoolArray)==1) + 1]);

% Find the end indices for each of StartInds. 
EndInds = find(diff(BoolArray)==-1);
if (numel(StartInds) ~= numel(EndInds))
    % This can be reached for a series that ends with 'true' on the last
    % index of BoolArray.
    EndInds = [EndInds; numel(BoolArray)];
end


end