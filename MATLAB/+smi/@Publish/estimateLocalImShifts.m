function [Shift, IntShift, ImageROIs, ImageStats] = ...
    estimateLocalImShifts(Image1, Image2, SubROISize, ...
    CorrParams, ShiftParams)
%estimateLocalImShifts estimates local shifts between two images.
%
% INPUT:
%   Image1: The stack to which Image2 is compared to, i.e. 
%           Stack1 is the reference stack. (mxn float array)
%   Image2: The stack for which the offset relative to Image1 
%           is to be determined. (mxn float array)
%   SubROISize: The size of local regions in which the shift will be
%               computed, ideally evenly divides [m, n].
%               (Pixels)(2x1 array)(Default = size(Image1))
%   CorrParams: Structure of parameters passed to smi_stat.findOffset().
%   ShiftParams: Structure of parameters passed to smi_stat.shiftImage()
%
% OUTPUT:
%   Shift: Sub-pixel offset of Image2 relative to Image1. (NROIsx2 array)
%   IntShift: Integer shift of Image2 relative to Image1. (NROIsx2 float)
%   ImageROIs: ROIs of the regions corresponding to the pixel offsets.
%              (NROIsx4 array)([YStart, XStart, YEnd, XEnd])
%   ImageStats: Structure containing some stats about the ImageROIs.
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
ImageSize = size(Image1);
if (~exist('SubROISize', 'var') || isempty(SubROISize))
    SubROISize = ImageSize.';
end
if (~exist('CorrParams', 'var') || isempty(CorrParams))
    CorrParams = struct([]);
end
if (~exist('ShiftParams', 'var') || isempty(ShiftParams))
    ShiftParams = struct([]);
end

% Split the images up into the sub-ROIs.
[DividedImages1, ImageROIs] = ...
    smi_helpers.subdivideImage(Image1, SubROISize);
DividedImages2 = smi_helpers.subdivideImage(Image2, SubROISize);

% Loop through each ROI and compute the local shift.
NROIs = size(ImageROIs, 1);
Shift = zeros(NROIs, 2);
IntShift = zeros(NROIs, 2);
for nn = 1:NROIs
    [ShiftCurrent, IntShiftCurrent] = smi_stat.findOffsetIter(...
        DividedImages1{nn}, DividedImages2{nn}, [], [], ...
        CorrParams, ShiftParams);
    Shift(nn, :) = ShiftCurrent(1:2);
    IntShift(nn, :) = IntShiftCurrent(1:2);
end

% If needed, compute some imaging stats. which might be useful to return.
if (nargin > 3)
    ImageStats.Image1Sums = cellfun(@(Image) sum(Image(:)), DividedImages1);
    ImageStats.Image2Sums = cellfun(@(Image) sum(Image(:)), DividedImages2);
end


end