function [SMD, ScaledData] = genLocalizations(obj)
%genLocalizations generates localizations from raw data.
% This method will generate localizations from an array of raw data.  This
% is done by first finding candidate ROIs in the raw data that may contain
% emitters, fitting a model function to the pixel values in that ROI, and
% then thresholding the resulting localizations.
% 
% REQUIRES:
%   DipImage, to use joinchannels()

% Created by:
%   David J. Schodt (Lidke Lab, 2020)
%       based on SMA_Core.fitStack written by Keith Lidke and 
%       Hanieh Mazloom-Farsibaf


% Construct an SMF structure from the class properties.
SMF = smi_core.SingleMoleculeFitting.createSMF();
SMF.Data.CameraType = obj.CameraType;
SMF.Data.CameraGain = obj.CameraGain;
SMF.Data.CameraOffset = obj.CameraOffset;
SMF.Data.CameraReadNoise = obj.CameraReadNoise;
SMF.FindROI.BoxSize = obj.BoxSize;
SMF.FindROI.PSFSigma = obj.PSFSigma;
SMF.FindROI.MinPhotons = obj.MinPhotons;

% Perform the gain and offset correction on the raw data.
[ScaledData] = smi_core.DataToPhotons.convertToPhotons(obj.RawData, ...
    obj.CameraGain, obj.CameraOffset, obj.CameraReadNoise);
                
% Generate candidate ROIs from the gain and offset corrected data.
FindROI = smi_core.FindROI(SMF, ScaledData);
[ROIStack, SMD] = FindROI.findROI();

% Pass the candidate ROIs to the fitting algorithm.  The output SMD from
% GaussMLE will contain localizations w.r.t. the ROIStack coordinates and
% thus we need to convert back to the full field of view before proceeding.
GaussMLE = smi_core.GaussMLE(SMF, ROIStack);
[SMD] = GaussMLE.gaussMLE(SMD);
SMD.X = SMD.X + SMD.XBoxCorner;
SMD.Y = SMD.Y + SMD.YBoxCorner;

% Set (but don't apply) the ThreshFlag in the SMD structure.
MinMax.X_SE = [0, SMF.Thresholding.MaxXY_SE];
MinMax.Y_SE = [0, SMF.Thresholding.MaxXY_SE];
if ~isempty(SMD.Z_SE)
    MinMax.Z_SE = [0, SMF.Thresholding.MaxZ_SE];
end
MinMax.PValue = [SMF.Thresholding.MinPValue, 1];
MinMax.PSFSigma = [SMF.Thresholding.MinPSFSigma, ...
    SMF.Thresholding.MaxPSFSigma];
MinMax.Photons = [SMF.Thresholding.MinPhotons, Inf];
MinMax.Bg = [0, SMF.Thresholding.MaxBG];
Threshold = smi_core.Threshold;
[SMD] = Threshold.setThreshFlag(SMD, MinMax);

% Store SMD and ScaledData as class properties.
obj.SMD = SMD;
obj.ScaledData = ScaledData;


end