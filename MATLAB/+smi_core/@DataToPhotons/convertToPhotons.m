function [Data, ReadNoise] = convertToPhotons(RawData, SMF, ...
    RawDataROI, CalibrationROI)
%convertToPhotons converts RawData to units of photons.
% This method will convert the image(s) contained in the array RawData to
% units of photons using the detector calibration information contained in
% the input CalibrationStruct.
%
% NOTE: The input RawData is assumed to be of the form 
%       RawData = Gain*RawData_photons + Offset where RawData_photons
%       is the output variable 'RawData' and RawData is given in units of
%       ADU.  It further assumed that one photoelectron correspond to one
%       photon, i.e., RawData_photons is really in units of photoelectrons.
%
% INPUTS:
%   RawData: An array whose entries correspond to pixels of the collection
%            camera. (ADU)(mxnxN numeric array)
%   SMF: A Single Molecule Fitting structure containing the following 
%        fields in a sub-structure 'Data' (see SingleMoleculeFitting class
%        for more details):
%        CameraGain: An array specifying the camera gain. This can be
%                    either a scalar (e.g., for a CCD) or a matrix (e.g., 
%                    for an sCMOS). 
%                    (ADU/e-)(1x1 array, or (m+a)x(n+b) array w/ a,b >= 0)
%        CameraOffset: An array specifying the camera offset.
%                      (ADU)(dimensions must match CameraGain)
%        CameraReadNoise: An array specifying the camera noise.
%                         (ADU^2)(dimensions must match CameraGain)
%   RawDataROI: The region of interest (ROI) of the input RawData on the
%               camera. 
%               (Pixels)([YStart, XStart, YEnd, XEnd])
%               (Default = [], meaning RawData is centered w.r.t. 
%               CameraGain/OffsetImage, or that CameraGain/OffsetImage are
%               scalars)
%   CalibrationROI: ROI of the CameraGain, OffsetImage, and CameraNoise (if
%                   they aren't scalars).  If CameraGain/OffsetImage aren't
%                   scalars, CalibrationROI must specify a region that at
%                   least contains the RawDataROI.
%                   (Pixels)([YStart, XStart, YEnd, XEnd])
%                   (Default = [1, 1, size(CameraGain)])
%
% OUTPUTS:
%   Data: The input RawData corrected for gain and offset (i.e., converted
%         to units of photons)(photons)(mxnxN numeric array)
%

% Created by:
%   David J. Schodt (Lidke Lab, 2020)
%       based on previous version by Hanieh Mazloom-Farsibaf


% Copy parameters from the CalibrationStruct to local variables and set
% defaults where appropriate.
if (~exist('RawDataROI', 'var') || isempty(RawDataROI))
    RawDataROI = [];
end
if (~exist('CalibrationROI', 'var') || isempty(CalibrationROI))
    CalibrationROI = [1, 1, size(SMF.Data.CameraGain)];
end

% Typecast arrays where appropriate.
RawData = single(RawData);
CameraGain = single(SMF.Data.CameraGain);
CameraOffset = single(SMF.Data.CameraOffset);
CameraReadNoise = single(SMF.Data.CameraReadNoise);

% Compare array sizes and ROI definitions to ensure consistency.
SizeRawData = size(RawData);
SizeCameraGain = size(CameraGain);
SizeCameraOffset = size(CameraOffset);
SizeCameraNoise = size(CameraReadNoise);
assert(all(SizeCameraGain == SizeCameraOffset), ...
    'DataToPhotons.convertToPhotons(): Input ''%s'' fields ', ...
    '''GainImage'' and ''OffsetImage'' must be the same size', ...
    inputname(2))
assert(...
    (prod(SizeCameraGain)==1) || all(SizeRawData(1:2)<=SizeCameraGain), ...
    ['DataToPhotons.convertToPhotons(): Input ''%s'' fields ', ...
    '''GainImage'' and ''OffsetImage'' must be either scalars or ', ...
    'they must be smaller than images in the image stack ''%s''.'], ...
    inputname(2), inputname(1))
assert(isempty(RawDataROI) ...
    || all((RawDataROI(3:4)-RawDataROI(1:2)+1)==SizeRawData(1:2)), ...
    ['DataToPhotons.convertToPhotons(): RawDataROI is not consistent ', ...
    'with the size of RawData'])
assert(isempty(RawDataROI) ...
    || all((RawDataROI(1:2)>=CalibrationROI(1:2)) ...
    & (RawDataROI(3:4)<=CalibrationROI(3:4))), ...
    ['DataToPhotons.convertToPhotons(): Input CalibrationROI', ...
    ' must contain the RawDataROI.'])
assert(...
    isempty(CameraReadNoise) || all(SizeCameraNoise == SizeCameraGain), ...
    ['DataToPhotons.convertToPhotons(): Input ''%s'' field ', ...
    '''CameraNoise'' must be either empty or the same size as ''%s'' ', ...
    ' field ''CameraGain'''], inputname(2), inputname(2))

% Isolate the portions of GainImage and OffsetImage corresponding to the
% specified RawDataROI of RawData (if needed).
if (prod(SizeCameraGain) > 1)
    % If RawDataROI is empty, we will assume RawData corresponds to the
    % center ROI of GainImage and OffsetImage.
    if isempty(RawDataROI)
        RawDataROI = repmat(CalibrationROI(1:2), [1, 2]) ...
            + [round((SizeCameraGain(1:2)-SizeRawData(1:2))/2), ...
            SizeRawData(1:2) + 1];
    end
    
    % Isolate the desired portions of GainImage and OffsetImage.
    RowIndices = RawDataROI(1):RawDataROI(3);
    ColumnIndices = RawDataROI(2):RawDataROI(4);
    CameraGainSub = CameraGain(RowIndices, ColumnIndices);
    CameraOffsetSub = CameraOffset(RowIndices, ColumnIndices);
    CameraNoiseSub = CameraReadNoise(RowIndices, ColumnIndices);
else
    CameraGainSub = CameraGain;
    CameraOffsetSub = CameraOffset;
    CameraNoiseSub = CameraReadNoise;
end

% Perform the gain and offset correction to convert RawData to units of
% photons.
% NOTE: MATLAB will still do this element-wise calculation correctly even
%       when size(RawData, 3) > 1 .
Data = (RawData-CameraOffsetSub) ./ CameraGainSub;
ReadNoise = CameraNoiseSub ./ (CameraGainSub.^2);

% Remove negative values from Data.
% NOTE: The old version of this code used a lower threshold of 0.01 instead
%       of 0.  I'm just using 0 because I don't see the benefit of using
%       0.01 (although one might exist that I'm just not seeing!).
Data(Data < 0) = 0; 


end