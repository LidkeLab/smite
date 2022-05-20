function [ROIStack, CameraGain, CameraOffset, CameraReadNoise] = ...
    extractROIs(SMD, SMF, CorrectData)
%extractROIs reloads raw data and extracts ROIs corresponding to SMD locs.
% This method reloads raw data defined by SMF.Data.FileDir and
% SMF.Data.FileName{1} and cuts out the ROIs defined by SMD.XBoxCorner,
% SMD.YBoxCorner, and SMF.BoxFinding.BoxSize.

% INPUTS:
%   SMD: Single Molecule Data Structure (see smi_core.SingleMoleculeData)
%   SMF: Single Molecule Fitting Structure (see
%        smi_core.SingleMoleculeFitting)
%   CorrectData: Flag indicating we should apply gain/offset correction to
%                the raw data before extracting the ROIs. (Default = true)
%
% OUTPUTS:
%   ROIStack: Stack of sub-ROIs extracted from the raw data.  The third
%             dimension shares the same indexing as SMD.XBoxCorner.
%             (SMF.BoxFinding.BoxSize x SMF.BoxFinding.BoxSize
%              x length(SMD.XBoxCorner) numeric array)
%   CameraGain/Offset/ReadNoise: Returned as a convenience since they're
%                                used internally (see
%                                smi_core.DataToPhotons).  Note that
%                                CameraReadNoise is the UN-CORRECTED
%                                readnoise image (i.e., it has units of
%                                ADU^2)

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Set defaults if needed.
if (~exist('CorrectData', 'var') || isempty(CorrectData))
    CorrectData = true;
end

% Load the calibration data.
if CorrectData
    DTP = smi_core.DataToPhotons(SMF);
else
    DTP = smi_core.DataToPhotons();
    DTP.CameraGain = 1;
    DTP.CameraOffset = 0;
end

% Loop over datasets in SMD, load data, and extract ROIs.
BoxSize = SMF.BoxFinding.BoxSize;
ROIStack = NaN(BoxSize, BoxSize, numel(SMD.XBoxCorner));
LD = smi_core.LoadData;
ROIStackInd = 1; % index to keep track of where to add new data in ROIStack
for nn = unique(SMD.DatasetNum.')
    % Load the raw data.
    [~, RawData] = LD.loadRawData(SMF, nn);

    % Apply gain/offset correction (note settings of DTP.CameraGain = 1 and
    % DTP.CameraOffset = 0 when CorrectData = false).
    DTP.RawData = RawData;
    CorrectedData = DTP.convertData();

    % Extract the ROIs.
    % NOTE: This code relies on SMD.XBoxCorner and SMD.YBoxCorner being
    %       defined such that the box length defined by
    %       SMF.BoxFinding.BoxSize does not push us beyond the border of
    %       the data (e.g., if the data is 256x256 and BoxSize=7,
    %       max(SMD.XBoxCorner)<=250 and max(SMD.YBoxCorner)<=250).
    SMDSub = ...
        smi_core.SingleMoleculeData.isolateSubSMD(SMD, SMD.DatasetNum==nn);
    for rr = 1:numel(SMDSub.XBoxCorner)
        ROIStack(:, :, ROIStackInd) = CorrectedData(...
            SMDSub.YBoxCorner(rr):(SMDSub.YBoxCorner(rr)+BoxSize-1), ...
            SMDSub.XBoxCorner(rr):(SMDSub.XBoxCorner(rr)+BoxSize-1), ...
            SMDSub.FrameNum(rr));
        ROIStackInd = ROIStackInd + 1;
    end
end

% Return extra outputs if needed.
if (nargout > 1)
    CameraGain = DTP.CameraGain;
    CameraOffset = DTP.CameraOffset;
    CameraReadNoise = DTP.CameraReadNoise;
end


end