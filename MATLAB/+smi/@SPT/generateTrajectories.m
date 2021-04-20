function generateTrajectories(obj)
%generateTrajectories generates trajectories from localizations in obj.SMD
% This method will generate trajectories from the localizations in obj.SMD
% using the parameters present in obj.SMF.

% Created by:
%   Hanieh Mazloom-Farsibaf (Lidke Lab, unknown date) as track.m
%   Reorganized and added comments, David J. Schodt (Lidke Lab, 2020)
%   Revised for smite, David J. Schodt (Lidke Lab, 2021)


% Set obj.DataROI if needed (this is not explicitly needed for tracking,
% but it's nice to keep it around in saved results).
if isempty(obj.DataROI)
    warning(['smi.SPT.generateTrajectories(): obj.DataROI property ', ...
        'not set. Default set to obj.DataROI = ', ...
        '[1, 1, obj.SMD.YSize, obj.SMD.XSize]'])
    obj.DataROI = [1, 1, obj.SMD.YSize, obj.SMD.XSize];
end
obj.SMD.DataROI = obj.DataROI;

% Estimate the density of off emitters (if requested).
if obj.EstimateRhoFromData
    RhoOnMean = mean(smi_core.SingleMoleculeData.computeDensity(obj.SMD));
    obj.SMF.Tracking.Rho_off = RhoOnMean ...
        * (obj.SMF.Tracking.K_off/obj.SMF.Tracking.K_on);
end

% Perform the frame-to-frame connection of the localizations.
for ff = min(obj.SMD.FrameNum):(max(obj.SMD.FrameNum)-1)
    % Create the frame-to-frame connection cost matrix.
    CostMatrix = smi.SPT.createCostMatrixFF(obj.SMD, obj.SMF, ...
        ff, obj.NonlinkMarker);

    % Perform the linear assignment problem to determine how we should link
    % together trajectories.
    Link12 = obj.solveLAP(CostMatrix);
    obj.SMD = obj.connectTrajFF(obj.SMD, Link12, ff);
end

% Perform the gap closing on the trajectory segments.
CostMatrix = obj.createCostMatrixGC(obj.SMD, obj.SMF, ...
    obj.NonlinkMarker, obj.UseSparseMatrices);
Link12 = obj.solveLAP(CostMatrix);
obj.SMD = obj.connectTrajGC(obj.SMD, Link12);

% Add the framerate and pixel size to the outputs.
% NOTE: This is a short term solution for back compatability.  WE WANT TO
%       REMOVE THIS EVENTUALLY!
[obj.SMD.FrameRate] = deal(obj.SMF.Data.FrameRate);
[obj.SMD.PixelSize] = deal(obj.SMF.Data.PixelSize);
obj.TR = smi_core.TrackingResults.convertSMDToTR(obj.SMD);


end