function [Mask, Shifts, ImageROIs] = ...
    genBFMask(FocusImageStruct, RefImage, MaxBrightfieldShift)
%genMaskedOverlays masks overlays based on brightfield shifts.
% This method loads brightfield images from the data files, computes shifts
% between those images, and then defines a mask based on those shifts.

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Loop through the provided FocusImageStructs and generate the masks.
if (isstruct(FocusImageStruct) && ~isempty(RefImage))
    % In this usage, each focus image is compared to the reference (as
    % opposed to focus images from another label).
    if isinf(MaxBrightfieldShift)
        Mask = ones(size(RefImage), 'logical');
        Shifts = NaN;
        ImageROIs = [1, 1, size(RefImage, 1:2)];
        return
    end
    [Shifts, ~, ImageROIs] = smi.Publish.computeBFShifts(...
        FocusImageStruct, RefImage);

    % Define the brightfield mask.
    Mask = smi.Publish.defineShiftMask(median(Shifts, 2), ...
        ImageROIs, ...
        MaxBrightfieldShift);
elseif iscell(FocusImageStruct)
    % In this usage, 'FocusImageStructs' is actually a cell array of
    % structs corresponding to two labels which will be compared to one
    % another.
    ImSize = size(FocusImageStruct{end}(1).Data.PreSeqImages, 1:2);
    Mask = ones([ImSize, size(FocusImageStruct, 1)], 'logical');
    Shifts = [];
    if isinf(MaxBrightfieldShift)
        ImageROIs = [1, 1, size(RefImage, 1:2)]; % define ImageROIs (MJW)
        return
    end
    for ii = 1:size(FocusImageStruct, 1)
        % Compute the local brightfield shifts from label 1 to 2.
        if (isempty(FocusImageStruct{ii, 1}) ...
                || isempty(FocusImageStruct{ii, 2}))
            continue
        end
        [LocalShiftMag, ~, ImageROIs] = smi.Publish.computeBFShifts(...
            FocusImageStruct{ii, 1}, FocusImageStruct{ii, 2});
        Shifts = cat(3, Shifts, LocalShiftMag);
        LocalShiftMagMedian = median(LocalShiftMag, 2);

        % Define the brightfield mask.
        Mask(:, :, ii) = smi.Publish.defineShiftMask(LocalShiftMagMedian, ...
            ImageROIs, ...
            MaxBrightfieldShift);
    end
end


end
