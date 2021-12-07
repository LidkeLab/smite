function [SMD, BestRegInd] = shiftToBestReg(SMD, RefImage, FocusImages)
%shiftToBestReg shifts coordinates based on best alignment results.
% This method shifts the coordinates SMD.X and SMD.Y based on the "best"
% results in brightfield registration as determined from
% AlignResultsStruct.
%
% INPUTS:
%   SMD: Single Molecule Data structure of localizations.
%   RefImage: Reference image at the focus. (YxX)
%   FocusImages: Structure containing the YxX images at the focal plane
%                taken before each SR sequence. (NDatasetsx1 struct)
%
% OUTPUTS:
%   SMD: Input SMD with SMD.X and SMD.Y shifted based on info. in
%        AlignResultsStruct.
%   BestRegInd: Dataset index selected as the "best" registration.

% Created by:
%   David J. Schodt (Lidke Lab 2021)


% Determine which dataset had the best channel registration (i.e., minimum
% XY shift in the focus image).
Shift = zeros(2, numel(FocusImages));
Params.SuppressWarnings = true;
for ii = 1:numel(FocusImages)
    CurrentShift = smi_stat.findOffsetIter(...
        median(FocusImages(ii).Data.PreSeqImages, 3), ...
        RefImage, [], [], Params);
    Shift(:, ii) = CurrentShift(1:2);
end
[~, BestRegInd] = min(sum(Shift.^2, 1));

% Make the drift correction reference to dataset BestRegInd.
SMD = smi_core.DriftCorrection.changeInterRef(SMD, BestRegInd);

% Shift each dataset based on the shift estimated from brightfield.
for nn = 1:SMD.NDatasets
    CurrentDSBool = (SMD.DatasetNum == nn);
    SMD.Y(CurrentDSBool) = SMD.Y(CurrentDSBool) - Shift(1, nn);
    SMD.X(CurrentDSBool) = SMD.X(CurrentDSBool) - Shift(2, nn);
end


end