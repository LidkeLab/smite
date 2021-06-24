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

% Loop through datasets and shift coordinates to the new reference.
% NOTE: The strange indexing used was done to improve speed.
NLocsPerDataset = groupcounts(SMD.DatasetNum);
NLocsCumulativeDS = [0; cumsum(NLocsPerDataset)];
XShift = SMD.DriftX(1, RefDatasetNum);
YShift = SMD.DriftY(1, RefDatasetNum);
for nn = 1:SMD.NDatasets
    CurrentDSInd = (1:NLocsPerDataset(nn)) + NLocsCumulativeDS(nn);
    [NLocsPerFrame, Frames] = groupcounts(SMD.FrameNum(CurrentDSInd));
    NLocsCumulativeFN = [0; cumsum(NLocsPerFrame)];
    for ff = 1:numel(Frames)
        CurrentFrameInd = (1:NLocsPerFrame(ff)) + NLocsCumulativeFN(ff);
        SMD.X(CurrentFrameInd) = SMD.X(CurrentFrameInd) + XShift;
        SMD.Y(CurrentFrameInd) = SMD.Y(CurrentFrameInd) + YShift;
    end
    SMD.DriftX(:, nn) = SMD.DriftX(:, nn) - XShift;
    SMD.DriftY(:, nn) = SMD.DriftY(:, nn) - XShift;
end


end