function [TR, SMD] = performFullAnalysis(obj)
%performFullAnalysis fits and tracks data pointed to by obj.SMF
% This method is the main run method for the smi.SPT class, meaning that it
% will load raw data, perform gain/offset correction, fit localizations to
% the data, create trajectories from the localizations, and then save the
% results.
%
% OUTPUTS:
%   TR: Tracking Results structure (see smi_core.TrackingResults)
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Prepare an SMLM class (we'll use this to load data, perform gain/offset
% correction, and fit the data).
if obj.SMF.Tracking.TryLowPValueLocs
    obj.SMFCopy = smi_core.SingleMoleculeFitting.reloadSMF(obj.SMF);
    obj.SMF.Thresholding.MinPValue = 0;
end
SMLM = smi.SMLM(obj.SMF);
SMLM.Verbose = obj.Verbose;

% Load data, perform gain/offset correction, and fit the data.
SMLM.analyzeAll()
obj.SMD = SMLM.SMD;
obj.SMDPreThresh = SMLM.SMDPreThresh;

% Perform the tracking.
obj.autoTrack()

% If a channel registration file is provided, apply the transform to our
% tracking results.
if ~isempty(obj.SMF.Data.RegistrationFilePath)
    if isfile(obj.SMF.Data.RegistrationFilePath)
        % Load the registration transform and apply it to obj.SMD and
        % recreate obj.TR based on the transformed SMD.
        load(obj.SMF.Data.RegistrationFilePath, 'RegistrationTransform')
        obj.SMDPreThresh = smi_core.ChannelRegistration.transformSMD(...
            RegistrationTransform, obj.SMDPreThresh);
        obj.SMD = smi_core.ChannelRegistration.transformSMD(...
            RegistrationTransform, obj.SMD);
        obj.TR = smi_core.TrackingResults.convertSMDToTR(obj.SMD);
    elseif (obj.Verbose > 0)
        warning(['smi.SPT.generateTrajectories(): the specified ', ...
            'SMF.Data.RegistrationFilePath cannot be found.'])
    end
end

% Remove short trajectories from the TR structure.
% NOTE: I'm leaving everything in SMD.  It might be nice to also threshold
%       short trajectories in SMD, but for now I'll leave it this way.
obj.TR = smi_core.TrackingResults.threshTrajLength(obj.TR, ...
    obj.SMF.Tracking.MinTrackLength);

% Make copies of TR and SMD for the outputs.
if nargout
    TR = obj.TR;
    SMD = obj.SMD;
end

% Save the results.
obj.saveResults()


end