function generateTrajectories(obj)
%generateTrajectories generates trajectories from localizations in obj.SMD
% This method will generate trajectories from the localizations in obj.SMD
% using the parameters present in obj.SMF.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Perform the frame-to-frame connection of the localizations.
obj.SMD = obj.genTrajFF(obj.SMD, obj.SMF, obj.RhoOff, obj.NonLinkMarker);

% Throw away unconnected low p-value localizations.
if obj.SMF.Tracking.TryLowPValueLocs
    obj.SMF.Thresholding.MinPValue = obj.SMFCopy.Thresholding.MinPValue;
    [GroupCounts, GroupIDs] = groupcounts(obj.SMD.ConnectID);
    UniqueTraj = GroupIDs(GroupCounts == 1);
    obj.SMD.ThreshFlag = ...
        int32(((obj.SMD.PValue<obj.SMF.Thresholding.MinPValue) ...
        & ismember(obj.SMD.ConnectID, UniqueTraj)));
    obj.SMD = smi_core.Threshold.applyThresh(obj.SMD, obj.Verbose);
end
obj.SMD.ConnectID = smi_helpers.compressToRange(obj.SMD.ConnectID);

% Perform the gap closing on the trajectory segments.
obj.SMD = obj.genTrajGC(obj.SMD, obj.SMF, obj.RhoOff, ...
    obj.NonLinkMarker, obj.UseSparseMatrices);

% Convert obj.SMD to a TR structure.
obj.TR = smi_core.TrackingResults.convertSMDToTR(obj.SMD);


end