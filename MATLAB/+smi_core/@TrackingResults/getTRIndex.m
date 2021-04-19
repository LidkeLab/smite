function [TRIndex] = getTRIndex(TR, TrajectoryIDs)
%getTRIndex returns the indices in TR corresponding to TrajectoryIDs.
% This method determines the indices in TR which contain info. about the
% trajectories 'TrajectoryIDs'.  For example, If TR(n).TrajectoryID = m, 
% getTRIndex(TR, m) will return n.
% 
% INPUTS:
%   TR: Tracking Results structure.
%   TrajectoryIDs: Array of trajectory IDs present in TR which will be
%                  joined together. (integer array)
%
% OUTPUTS:
%   TRIndex: Indices in the input TR structure corresponding to the
%            input TrajectoryIDs. If the requested TrajectoryID wasn't
%            found, the corresponding index is returned as NaN. 
%            (integer array, same size as TrajectoryIDs)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Ensure TrajectoryIDs is a row vector (see usage in find() below).
if iscolumn(TrajectoryIDs)
    TrajectoryIDs = TrajectoryIDs.';
end

% Find the requested indices in TR.
% NOTE: We need the ~ for the second output of find() because requesting
%       that output produces the desired behavior.  If it's removed, find()
%       returns a "column stacked" index which we'd have to convert
%       ourselves.
TrajectoryIDAll = cell2mat({TR.TrajectoryID}.');
[TRIndexTemp, ~] = find(TrajectoryIDAll == TrajectoryIDs);
TRIndex = NaN(numel(TrajectoryIDs), 1);
TRIndex(ismember(TrajectoryIDs, TrajectoryIDAll)) = TRIndexTemp;


end