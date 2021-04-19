function [TR] = joinTraj(TR, TrajectoryIDs)
%joinTraj joins a set of trajectories into one trajectory.
% This method joins the trajectorys with IDs 'TrajectoryIDs' into a single
% trajectory with the new ID being min(TrajectoryIDs).
%
% INPUTS:
%   TR: Tracking Results structure.
%   TrajectoryIDs: Array of trajectory IDs present in TR which will be
%                  joined together. (integer array)
%
% OUTPUTS:
%   TR: Tracking Results Structure with specified trajectories combined
%       into one trajectory.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Remove duplicates and sort the TrajectoryIDs input.
TrajectoryIDs = sort(unique(TrajectoryIDs), 'ascend');

% Loop through each of TrajectoryIDs and begin merging together.
TRIndices = smi_core.TrackingResults.getTRIndex(TR, TrajectoryIDs);
for ii = 2:numel(TrajectoryIDs)
    if any(intersect(TR(TRIndices(1)).FrameNum, ...
            TR(TRIndices(ii)).FrameNum))
        warning(['Trajectories %i and %i cannot be joined because ', ...
            'they were observed in the same frame'], ...
            TR(TRIndices(1)).TrajectoryID, TR(TRIndices(ii)).TrajectoryID)
        continue
    end
    TR(TRIndices(1)).X = [TR(TRIndices(1)).X; TR(TRIndices(ii)).X];
    TR(TRIndices(1)).Y = [TR(TRIndices(1)).Y; TR(TRIndices(ii)).Y];
    TR(TRIndices(1)).X_SE = [TR(TRIndices(1)).X_SE; ...
        TR(TRIndices(ii)).X_SE];
    TR(TRIndices(1)).Y_SE = [TR(TRIndices(1)).Y_SE; ...
        TR(TRIndices(ii)).Y_SE];
    TR(TRIndices(1)).FrameNum = [TR(TRIndices(1)).FrameNum; ...
        TR(TRIndices(ii)).FrameNum];
    TR(TRIndices(1)).Photons = [TR(TRIndices(1)).Photons; ...
        TR(TRIndices(ii)).Photons];
    TR(TRIndices(1)).Bg = [TR(TRIndices(1)).Bg; ...
        TR(TRIndices(ii)).Bg];
    TR(TRIndices(1)).IndSMD = [TR(TRIndices(1)).IndSMD; ...
        TR(TRIndices(ii)).IndSMD];
    TR(TRIndices(ii)) = [];
end

% Sort the resulting trajectory w.r.t. FrameNum.
[TR(TRIndices(1)).FrameNum, SortIndices] = sort(TR(TRIndices(1)).FrameNum);
TR(TRIndices(1)).X = TR(TRIndices(1)).X(SortIndices);
TR(TRIndices(1)).Y = TR(TRIndices(1)).Y(SortIndices);
TR(TRIndices(1)).X_SE = TR(TRIndices(1)).X_SE(SortIndices);
TR(TRIndices(1)).Y_SE = TR(TRIndices(1)).Y_SE(SortIndices);
TR(TRIndices(1)).Photons = TR(TRIndices(1)).Photons(SortIndices);
TR(TRIndices(1)).Bg = TR(TRIndices(1)).Bg(SortIndices);
TR(TRIndices(1)).IndSMD = TR(TRIndices(1)).IndSMD(SortIndices);



end