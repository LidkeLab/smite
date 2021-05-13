function [KOn, KOff, KBleach] = estimateRateParameters(SMD, ...
    MinRate, MaxRate)
%estimateRateParameters estimates blinking rates for localizations in SMD
% This method will make an estimate of the blinking kinetics (KOn, KOff,
% and KBleach) based on the (pre-clustered) localizations in SMD.
%
% INPUTS:
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)
%        containing localizations associated by SMD.ConnectID
%
% OUTPUTS:
%   KOn: Rate parameter for blinking on (between MinRate and MaxRate)
%        (1 / frame)
%   KOff: Rate parameter for blinking off 
%         (between MinRate and MaxRate)(1 / frame)
%   KBleach: Rate parameter for photobleaching
%            (between MinRate and MaxRate)(1 / frame)
%   MinRate: Minimum rate parameter allowed (Default = 1e-5)
%   MaxRate: Maximum rate parameter allowed (Default = 1-MinRate)
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('MinRate', 'var') || isempty(MinRate))
    MinRate = 1e-5;
end
if (~exist('MaxRate', 'var') || isempty(MaxRate))
    MaxRate = 1 - MinRate;
end

% Grab some arrays from the SMD (sometimes can help with speed).
FrameNum = SMD.FrameNum;
NFrames = SMD.NFrames;
ConnectID = SMD.ConnectID;

% Determine the bright durations, dark durations, and total lifetime of
% each observed cluster (assuming it's one emitter blinking on and off and
% then bleaching).
OnDurations = [];
OffDurations = [];
OnDurationsTotal = [];
NClusters = numel(unique(ConnectID));
for nn = 1:NClusters
    % NOTE: I'm taking unique(FrameNum) here because our pre-clustering
    %       allows multiple localizations per frame.
    CurrentClusterBool = (ConnectID == nn);
    FrameDiff = diff(unique(FrameNum(CurrentClusterBool)));
    OffIndices = find(FrameDiff ~= 1);
    NOffEvents = numel(OffIndices);
    if (NOffEvents < 2)
        % If there was only 1 off event, we should play it safe and assume
        % it was a bleaching event.
        continue
    end

    % For the on durations, we want to skip the last duration since that
    % might be a bleaching event. We'll estimate the total on duration by
    % counting the total number of localizations in the cluster (since each
    % one is in a distinct frame, thus that's roughly the number of frames
    % the emitter was on) unless the emitter was observed in the last
    % frame.
    OffDurations = [OffDurations; FrameDiff(OffIndices) - 1];
    OnDurations = [OnDurations; ...
        OffIndices(1); diff(OffIndices(1:NOffEvents-1))];
    if (NFrames > max(FrameNum(CurrentClusterBool)))
        OnDurationsTotal = [OnDurationsTotal; sum(CurrentClusterBool)];
    end
end
KOn = max(MinRate, min(MaxRate, -log(1 - 1/mean(OffDurations))));
if isnan(KOn)
    KOn = 1;
end
KOff = max(MinRate, min(MaxRate, -log(1 - 1/mean(OnDurations))));
if isnan(KOff)
    KOff = 1;
end
KBleach = max(MinRate, min(MaxRate, -log(1 - 1/mean(OnDurationsTotal))));
if isnan(KBleach)
    KBleach = 1;
end


end