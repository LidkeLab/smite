function [LocalShiftMag, LocalShift1To2, ImageROIs] = ...
    computeBFShifts(FocusImStruct1, FocusImStruct2, SubROISize)
%computeBFShifts computes subROI shifts from brightfield data.
% This method will compute sub-ROI shifts between the brightfield images
% present in FocusImStruct1 and 2.
%
% INPUT:
%   FocusImStruct1: Structure containing the focus images for the first
%                   label (see the loaded structure BFStruct in, e.g.,
%                   smi.Publish.processLabel()).
%   FocusImStruct2: Structure containing the focus images for the second
%                   label (see the loaded structure BFStruct in, e.g.,
%                   smi.Publish.processLabel()).
%   SubROISize: Size of the subregions. (pixels)(Default = [32; 32])
%   ImageROIs: ROIs of the regions corresponding to the pixel offsets.
%              (NROIsx4 array)([YStart, XStart, YEnd, XEnd])
%
% OUTPUT:
%   LocalShiftMag: Magnitude of the shift differences between labels 1 and
%                  2. (pixels)
%   LocalShift1To2: Vector shifts between labels 1 and 2.
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Set defaults.
if (~exist('SubROISize', 'var') || isempty(SubROISize))
    SubROISize = [32; 32];
end

% Loop through datasets and compute the shifts between the labels.
% NOTE: For now, I'm using the very last image taken just before the
%       sequence starts. It might be useful to use the median over the
%       PreSeqImages, but I'm avoiding that just in case there is a
%       dramatic sample drift.
CorrParams.SuppressWarnings = true;
LocalShift1To2 = ...
    zeros(prod(size(FocusImStruct1(1).Data.PreSeqImages(:, :, end)) ...
    ./ SubROISize.'), numel(FocusImStruct1), 2);
for ii = 1:numel(FocusImStruct1)
    [LocalShift1To2(:, ii, :), ~, ImageROIs] = ...
        smi.Publish.estimateLocalImShifts(...
        FocusImStruct1(ii).Data.PreSeqImages(:, :, end), ...
        FocusImStruct2(ii).Data.PreSeqImages(:, :, end), ...
        SubROISize, CorrParams);
end
LocalShiftMag = sqrt(sum(LocalShift1To2.^2, 3));


end