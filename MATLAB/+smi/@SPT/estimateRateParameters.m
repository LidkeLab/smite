function [KOn, KOff] = estimateRateParameters(SMD)
%estimateRateParameters estimates blinking rates for localizations in SMD
% This method will make an estimate of the blinking kinetics (KOn, KOff) 
% based on the localizations in SMD associated by SMD.ConnectID.
%
% INPUTS:
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)
%        containing localizations associated by SMD.ConnectID
%
% OUTPUTS:
%   KOn: Rate parameter for blinking on (1 / frame)
%   KOff: Rate parameter for blinking off (1 / frame)
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Grab some arrays from the SMD (sometimes can help with speed).
FrameNum = SMD.FrameNum;
ConnectID = SMD.ConnectID;

% Determine the bright and dark durations for each trajectory.
OnDurations = [];
OffDurations = [];
UniqueIDs = unique(ConnectID);
for ii = 1:numel(UniqueIDs)
    CurrentTraj = (SMD.ConnectID == UniqueIDs(ii));
    FrameDiff = diff(FrameNum(CurrentTraj));
    OffIndices = find(FrameDiff ~= 1);
    if isempty(OffIndices)
        continue
    end
    OffDurations = [OffDurations; FrameDiff(OffIndices) - 1];
    OnDurations = [OnDurations; OffIndices(1); diff(OffIndices)];
end

% Make an estimate of the rate parameters, setting a default value if not
% enough information is available.
KOn = -log(1 - 1/mean(OffDurations));
KOff = -log(1 - 1/mean(OnDurations));
if isnan(KOn)
    KOn = 1;
end
if isnan(KOff)
    KOff = 0;
end


end