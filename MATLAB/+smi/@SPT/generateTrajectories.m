function generateTrajectories(obj)
%generateTrajectories generates trajectories from localizations in obj.SMD
% This method will generate trajectories from the localizations in obj.SMD
% using the parameters present in obj.SMF.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Perform the frame-to-frame connection of the localizations.
obj.SMD.ConnectID = (1:numel(obj.SMD.FrameNum)).';
for ff = min(obj.SMD.FrameNum):(max(obj.SMD.FrameNum)-1)
    % Create the frame-to-frame connection cost matrix.
    CostMatrix = smi.SPT.createCostMatrixFF(obj.SMD, obj.SMF, ...
        obj.DiffusionConstant, ff, obj.NonlinkMarker);
    if (numel(CostMatrix) < 2)
        % If there's only one localization considered, there's no use in
        % proceeding.
        continue
    end
    
    % Perform the linear assignment problem to determine how we should link
    % together trajectories.
    Link12 = obj.solveLAP(CostMatrix);
    obj.SMD = obj.connectTrajFF(obj.SMD, Link12, ff);
end

% Throw away unconnected low p-value localizations.
if obj.SMF.Tracking.TryLowPValueLocs
    obj.SMF.Thresholding.MinPValue = obj.SMFCopy.Thresholding.MinPValue;
    [GroupCounts, GroupIDs] = groupcounts(obj.SMD.ConnectID);
    UniqueTraj = GroupIDs(GroupCounts == 1);
    obj.SMD.ThreshFlag = ...
        ((obj.SMD.PValue<obj.SMF.Thresholding.MinPValue) ...
        & ismember(obj.SMD.ConnectID, UniqueTraj));
    obj.SMD = smi_core.Threshold.applyThresh(obj.SMD, obj.Verbose);
end
obj.SMD.ConnectID = smi.SPT.validifyConnectID(obj.SMD.ConnectID);

% Perform the gap closing on the trajectory segments.
CostMatrix = obj.createCostMatrixGC(obj.SMD, obj.SMF, ...
    obj.DiffusionConstant, obj.NonlinkMarker, obj.UseSparseMatrices);
Link12 = obj.solveLAP(CostMatrix);
obj.SMD = obj.connectTrajGC(obj.SMD, Link12);

% Convert obj.SMD to a TR structure.
obj.TR = smi_core.TrackingResults.convertSMDToTR(obj.SMD);


end