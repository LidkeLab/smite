function [Mask] = defineShiftMask(LocalImShifts, ImageROIs, MaxShift, ImSize)
%generateShiftMask generates a mask from the array of shifts LocalImShifts.
%
% INPUT:
%   LocalImShifts: Array of shifts computed in estimateLocalImShifts().
%   ImageROIs: ROIs corresponding to LocalImShifts output by
%              estimateLocalImShifts().
%   MaxShift: Maximum shift allowed by output mask. 
%             (same units as LocalImShifts)
%   ImSize: Size of the output mask. (Default determined from
%           LocalImShifts).
%
% OUTPUT:
%   Mask: Logical array to be used to mask images. (ImSize logical array)
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Set defaults.
ImageROISize = max(ImageROIs(:, 3:4));
if (~exist('ImSize', 'var') || isempty(ImSize))
    ImSize = ImageROISize;
end
if (~exist('MaxShift', 'var') || isempty(MaxShift))
    MaxShift = inf;
end

% Generate the binary mask.
Mask = zeros(ImageROISize, 'logical');
LocalShiftMag = sqrt(sum(LocalImShifts.^2, 2));
for ii = 1:size(LocalImShifts, 1)
    Mask(ImageROIs(ii, 1):ImageROIs(ii, 3), ImageROIs(ii, 2):ImageROIs(ii, 4)) = ...
        (LocalShiftMag(ii) <= MaxShift);
end

% If needed, upsample the mask.
Mask = imresize(Mask, max(ImSize ./ ImageROISize));


end