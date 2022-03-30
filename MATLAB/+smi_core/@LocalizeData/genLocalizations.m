function [SMD, SMDPreThresh] = genLocalizations(obj)
%genLocalizations generates localizations from scaled data.
% This method will generate localizations from an array of data by first
% first finding candidate ROIs in the data that may contain emitters, 
% fitting a model function to the pixel values in that ROI, and then 
% thresholding the resulting localizations.
% 
% OUTPUTS: 
%   SMD: Single Molecule Data structure (see SingleMoleculeData class) with
%        only valid localizations (i.e., those that passed all thresholds).
%   SMDPreThresh: SMD structure with all found localizations, even those
%                 that did not pass the thresholds.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)
%       based on SMA_Core.fitStack written by Keith Lidke and 
%       Hanieh Mazloom-Farsibaf


% Generate candidate ROIs from the gain and offset corrected data.
if (obj.Verbose > 1)
    fprintf(['\tLocalizeData.genLocalizations(): ', ...
        'Finding candidate ROIs from the input data...\n'])
elseif (obj.Verbose > 0)
    fprintf(['\tLocalizeData.genLocalizations(): ', ...
        'Generating localizations from the input data...\n'])
end
FindROI = smi_core.FindROI(obj.SMF, obj.ScaledData);
FindROI.Verbose = obj.Verbose;
[ROIStack, SMDCandidates] = FindROI.findROI();
if (obj.Verbose > 2)
    fprintf(['\t\tLocalizeData.genLocalizations(): ', ...
        '%i candidate ROIs found.\n'], numel(SMDCandidates.FrameNum))
end

% Pass the candidate ROIs to the fitting algorithm.  The output SMD from
% GaussMLE will contain localizations w.r.t. the ROIStack coordinates and
% thus we need to convert back to the full field of view before proceeding.
if (obj.Verbose > 1)
    fprintf(['\tLocalizeData.genLocalizations(): ', ...
        'Fitting candidate ROIs with smi_core.GaussMLE...\n'])
end
GaussMLE = smi_core.GaussMLE(obj.SMF, ROIStack);
[SMDCandidates] = GaussMLE.gaussMLE(SMDCandidates);
SMDCandidates.X = SMDCandidates.X + SMDCandidates.XBoxCorner;
SMDCandidates.Y = SMDCandidates.Y + SMDCandidates.YBoxCorner;

% Threshold localizations.
MinMax = smi_core.Threshold.setMinMax(obj.SMF);
Threshold = smi_core.Threshold;
Threshold.Verbose = obj.Verbose;
[SMDPreThresh] = Threshold.setThreshFlag(SMDCandidates, MinMax);
if obj.SMF.Thresholding.On
   [SMD] = Threshold.applyThresh(SMDPreThresh, obj.Verbose);
else
   SMD = SMDPreThresh;
end
obj.SMDPreThresh = SMDPreThresh;
obj.SMD = SMD;

% If needed, produce some color overlay plots.
% NOTE: Extra cases kept in case of later changes.
%switch obj.Verbose
%    case 0
%    case 1
%    case 2
%    case 3
%        obj.colorOverlay();
%    otherwise
%        obj.colorOverlay();
%end
if obj.Verbose >= 3
   if isempty(obj.ResultsDir)
      obj.colorOverlay();
   else
      ResultsDir = obj.ResultsDir;
      FileBaseName = 'LocalizeDataRGBImage';
      RGBImage = obj.colorOverlay();
      smi_helpers.writeMPEG4(ResultsDir, FileBaseName, RGBImage);
   end
end

end
