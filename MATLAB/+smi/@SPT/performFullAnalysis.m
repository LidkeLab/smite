function [TR, SMD, SMDPreThresh] = performFullAnalysis(obj)
%performFullAnalysis fits and tracks data pointed to by obj.SMF
% This method is the main run method for the smi.SPT class, meaning that it
% will load raw data, perform gain/offset correction, fit localizations to
% the data, create trajectories from the localizations, and then save the
% results.
%
% OUTPUTS:
%   TR: Tracking Results structure (see smi_core.TrackingResults)
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)
%   SMDPreThresh: Output 'SMD' before thresholds are applied.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Load data, perform gain/offset correction, fit the data, and threshold
% the fits.
if obj.SMF.Tracking.TryLowPValueLocs
    obj.SMFCopy = smi_core.SingleMoleculeFitting.reloadSMF(obj.SMF);
    obj.SMF.Thresholding.MinPValue = 0;
end
LD = smi_core.LoadData;
[~, RawData] = LD.loadRawData(obj.SMF, 1);
DTP = smi_core.DataToPhotons(obj.SMF, RawData, [], [], obj.Verbose);
obj.ScaledData = DTP.convertData();
LD = smi_core.LocalizeData(obj.ScaledData, obj.SMF, obj.Verbose);
[obj.SMD, obj.SMDPreThresh] = LD.genLocalizations();
obj.SMD.NDatasets = 1;
obj.SMD.DatasetNum = ones(size(obj.SMD.FrameNum));

% Add pixel size and framerate to SMD (this is a temporary workaround for
% other codes, e.g., the diffusion estimator, and should be removed later).
obj.SMD.PixelSize = obj.SMF.Data.PixelSize;
obj.SMD.FrameRate = obj.SMF.Data.FrameRate;
obj.SMDPreThresh.PixelSize = obj.SMF.Data.PixelSize;
obj.SMDPreThresh.FrameRate = obj.SMF.Data.FrameRate;

% Check if localizations were generated.  If none were generated, issue a
% warning and do not proceed.
if isempty(obj.SMD.FrameNum)
    if (obj.Verbose > 0)
        warning(['smi.SPT.performFullAnalysis(): No localizations were ', ...
            'generated from the provided dataset.'])
    end
    if nargout
        TR = obj.TR;
        SMD = obj.SMD;
    end
    return
end

% Perform the tracking.
obj.autoTrack()

% If a channel registration file is provided, apply the transform to our
% tracking results.
if ~isempty(obj.SMF.Data.RegistrationFilePath)
    if isfile(obj.SMF.Data.RegistrationFilePath)
        % Load the registration transform and determine which transform
        % should be used based on the provided ROIs.
        load(obj.SMF.Data.RegistrationFilePath, ...
            'RegistrationTransform', 'FiducialROI', 'RegistrationErrorLOO')
        MatchedROI = find(all(FiducialROI == ...
            repmat(obj.SMF.Data.DataROI, [size(FiducialROI, 1), 1]), 2), ...
            1);
        
        % Apply the transform to the localizations. If no transform was
        % matched, issue a warning and proceed.
        if isempty(MatchedROI)
            warning(['smi.SPT.generateTrajectories(): None of the ', ...
                'available transforms match SMF.Data.DataROI!'])
        else
            obj.SMDPreCR = obj.SMD;
            obj.SMDPreThreshPreCR = obj.SMDPreThresh;
            obj.SMDPreThresh = smi_core.ChannelRegistration.transformSMD(...
                RegistrationTransform{MatchedROI}, obj.SMDPreThresh);
            obj.SMDPreThresh.RegError = RegistrationErrorLOO(MatchedROI);
            obj.SMD = smi_core.ChannelRegistration.transformSMD(...
                RegistrationTransform{MatchedROI}, obj.SMD);
            obj.SMD.RegError = RegistrationErrorLOO(MatchedROI);
        end
        obj.TRPreCR = smi_core.TrackingResults.convertSMDToTR(obj.SMDPreCR);
        obj.TR = smi_core.TrackingResults.convertSMDToTR(obj.SMD);
    elseif (obj.Verbose > 0)
        warning(['smi.SPT.generateTrajectories(): The specified ', ...
            'SMF.Data.RegistrationFilePath cannot be found!'])
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
    SMDPreThresh = obj.SMDPreThresh;
end

% Save the results.
if ~obj.IsTestRun
    obj.saveResults()
end

% Share an update in Command Window.
if (obj.Verbose > 0)
    fprintf('smi.SPT.performFullAnalysis(): Finished tracking file %s\n', ...
        obj.SMF.Data.FileName{1})
end


end