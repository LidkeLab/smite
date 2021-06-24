function [SMD] = changeInterRef(SMD, RefDatasetNum)
%changeInterRef shifts SMD coordinates to new dataset reference.
% This method will change the "reference" dataset number for inter-dataset
% drift results.  Specifically, this method can be applied AFTER
% inter-dataset drift correction was performed on 'SMD', with the intention
% being to change the asssumed reference dataset of dataset 1 to the
% dataset 'RefDatasetNum'.  This method assumes that 'SMD' arrays are
% sorted w.r.t. DatasetNum and FrameNum.
%
% INPUTS:
%   SMD: Single Molecule Data structure with populated fields DriftX and
%        DriftY and array fields sorted w.r.t. DatasetNum and FrameNum.
%   RefDatasetNum: Dataset number that will be the new reference.
%
% OUTPUTS:
%   SMD: Input SMD with coordinates shifted to the new reference dataset.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Check if RefDatasetNum is a valid choice.
if (RefDatasetNum > SMD.NDatasets)
    warning('Input ''RefDatasetNum'' can''t exceed ''SMD.NDatasets''')
    return
end

% Shift coordinates to the new reference and update the drift arrays.
XShift = SMD.DriftX(1, RefDatasetNum);
YShift = SMD.DriftY(1, RefDatasetNum);
SMD.X = SMD.X + XShift;
SMD.Y = SMD.Y + YShift;
SMD.DriftX = SMD.DriftX - XShift;
SMD.DriftY = SMD.DriftY - YShift;


end