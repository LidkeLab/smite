function [Mask] = genBFMask(FocusImageStructs, MaxBrightfieldShift)
%genMaskedOverlays masks overlays based on brightfield shifts.
% This method loads brightfield images from the data files, computes shifts
% between those images, and then defines a mask based on those shifts.

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Loop through the provided FocusImageStructs and generate the masks.
ImSize = size(FocusImageStructs{end}(1).Data.PreSeqImages, 1:2);
Mask = ones([ImSize, size(FocusImageStructs, 1)], 'logical');
if isinf(MaxBrightfieldShift)
    return
end
for ii = 1:size(FocusImageStructs, 1)
    % Compute the local brightfield shifts from label 1 to 2.
    if (isempty(FocusImageStructs{ii, 1}) ...
            || isempty(FocusImageStructs{ii, 2}))
        continue
    end
    [LocalShiftMag, ~, ImageROIs] = smi.Publish.computeBFShifts(...
        FocusImageStructs{ii, 1}, FocusImageStructs{ii, 2});
    LocalShiftMagMedian = median(LocalShiftMag, 2);

    % Define the brightfield mask.
    Mask(:, :, ii) = smi.Publish.defineShiftMask(LocalShiftMagMedian, ...
        ImageROIs, ...
        MaxBrightfieldShift);    
end


end