function [PixelOffsets, SubPixelOffsets, ImageROIs] = ...
    estimateLocalShifts(Image1, Image2, SubROISize, MaxOffset, UseGPU)
%estimateLocalShifts estimates local shifts between two images.
%
% INPUT:
%   Image1: The stack to which Image2 is compared to, i.e. 
%           Stack1 is the reference stack. (mxn float array)
%   Image2: The stack for which the offset relative to Image1 
%           is to be determined. (mxn float array)
%   SubROISize: The size of local regions in which the shift will be
%               computed, ideally evenly divides [m, n].
%               (Pixels)(2x1 array)(Default = size(Image1))
%   MaxOffset: Max offset for which the cross correlation is computed
%              between the two images. 
%              (Pixels)(Default = ceil(SubROISize / 4))
%
% OUTPUT:
%   PixelOffset: The integer pixel offset of Image2 relative to Image1,
%                determined based on the location of the peak of the xcorr
%                coefficient between the two images. (NROIsx2 array)
%   SubPixelOffset: The sub-pixel offset of Image2 relative to Image1, 
%                   approximated based on a 2nd order polynomial fit(s) to 
%                   the cross-correlation. (NROIsx2 float)
%   ImageROIs: ROIs of the regions corresponding to the pixel offsets.
%              (NROIsx4 array)([YStart, XStart, YEnd, XEnd])
%
% REQUIRES: 
%   matlab-instrument-control, to use findStackOffset()
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
ImageSize = size(Image1);
if (~exist('SubROISize', 'var') || isempty(SubROISize))
    SubROISize = ImageSize.';
end
if (~exist('MaxOffset', 'var') || isempty(MaxOffset))
    MaxOffset = ceil(SubROISize / 4);
end
if (~exist('UseGPU', 'var') || isempty(UseGPU))
    UseGPU = logical(gpuDeviceCount());
end

% Split the images up into the sub-ROIs.
[DividedImages1, ImageROIs] = ...
    smi_helpers.subdivideImage(Image1, SubROISize);
DividedImages2 = smi_helpers.subdivideImage(Image2, SubROISize);

% Loop through each ROI and compute the local shift.
NROIs = size(ImageROIs, 1);
PixelOffsets = zeros(NROIs, 2);
SubPixelOffsets = PixelOffsets;
for nn = 1:NROIs
    [Offset, SubOffset] = MIC_Reg3DTrans.findStackOffset(...
        DividedImages1{nn}, DividedImages2{nn}, ...
        [MaxOffset, 1], [], [], 0, UseGPU);
    PixelOffsets(nn, :) = Offset(1:2);
    SubPixelOffsets(nn, :) = SubOffset(1:2);
end


end