function [TR, SMD] = generateTrajectories(obj)
%generateTrajectories generates trajectories from localizations in obj.SMD
% This method will generate trajectories from the localizations in obj.SMD
% using the parameters present in obj.SMF.
%
% OUTPUTS:
%   TR: Tracking Results structure (see smi_core.TrackingResults)
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)

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

% Grab a bunch of parameters from the SMF structure (this might help make
% the code more readable).
K_on = obj.SMF.Params_Track.K_on;
K_off = obj.SMF.Params_Track.K_off;
D = obj.SMF.Params_Track.D;
MaxFrameGap = obj.SMF.Params_Track.MaxFrameGap;
MaxDistGC = obj.SMF.Params_Track.MaxDistGC;
MaxDistFF = obj.SMF.Params_Track.MaxDistFF;

% Calculate the mean density of emitters in the data.
fprintf('Estimating the density of emitters in the data...\n')
obj.TD.TrajectoryID = [];
RhoOnMean = mean(obj.calcDensity(obj.TD)); % mean density of 'on' emitters
RhoOffMean = (K_off/K_on) * RhoOnMean; % mean density of 'off' emitters

% If the framerate wasn't provided, error out in a useful way.
if isempty(obj.TD.FrameRate)
    error('SMA_SPT.track(): FrameRate not set inside of the TD structure')
end

% Perform the frame-to-frame connection of the localizations.
fprintf('Performing the frame-to-frame connection process...\n')
for ff = min(obj.TD.FrameNum):(max(obj.TD.FrameNum)-1)
    % Create the frame-to-frame connection cost matrix.
    [CM] = obj.createCM_FF(obj.TD, ff, K_on, K_off, RhoOffMean, D, ...
        MaxDistFF);

    % Perform the linear assignment problem to determine how we should link
    % together trajectories.
    [Link12] = obj.solveLAP(CM);
    [obj.TD] = obj.connectFrame(obj.TD, Link12, ff);
end

% Perform the gap closing on the trajectory segments.
fprintf('Performing the gap closing process...\n')
[CM] = obj.createCM_GC(obj.TD, K_on, K_off, RhoOffMean, ...
    D, MaxFrameGap, MaxDistGC);
[Link12] = obj.solveLAP(CM);
[obj.TD] = obj.connectGaps(obj.TD, Link12);


end