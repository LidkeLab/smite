function [PixelOffsets, SubPixelOffsets, BlockROIs] = ...
    estimateLocalShifts(Image1, Image2, BlockSize, UseGPU)
%estimateLocalShifts estimates local shifts between two images.
%
% INPUT:
%   Image1: The stack to which Image2 is compared to, i.e. 
%           Stack1 is the reference stack. (mxn float array)
%   Image2: The stack for which the offset relative to Image1 
%           is to be determined. (mxn float array)
%   BlockSize: The size of local regions in which the shift will be
%              computed, ideally evenly divides the [m, n].
%              (Pixels)(2x1 array)(Default = size(Image1))
%
% OUTPUT:
%   PixelOffset: The integer pixel offset of Image2 relative to Image1,
%                determined based on the location of the peak of the xcorr
%                coefficient between the two images. (NROIsx2 array)
%   SubPixelOffset: The sub-pixel offset of Image2 relative to Image1, 
%                   approximated based on a 2nd order polynomial fit(s) to 
%                   the cross-correlation. (NROIsx2 float)
%   BlockROIs: ROIs of the regions corresponding to the pixel offsets.
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
if (~exist('BlockSize', 'var') || isempty(BlockSize))
    BlockSize = ImageSize.';
end
if (~exist('UseGPU', 'var') || isempty(UseGPU))
    UseGPU = logical(gpuDeviceCount());
end

% Ensure arrays are properly shaped.
if isrow(BlockSize)
    BlockSize = BlockSize.';
end

% Define the ROIs of the local regions which we'll find the shifts within.
NDivisions = ceil(size(Image1).' ./ BlockSize);
YStart = repmat(1 + BlockSize(1)*(0:(NDivisions(1)-1)).', ...
    [NDivisions(2), 1]);
XStart = repelem(1 + BlockSize(2)*(0:(NDivisions(2)-1)).', ...
    NDivisions(1));
BlockROIs = [YStart, XStart, ...
    min(ImageSize(1), YStart+BlockSize(1)-1), ...
    min(ImageSize(2), XStart+BlockSize(2)-1)];

% Loop through each ROI and compute the local shift.
NROIs = prod(NDivisions);
PixelOffsets = zeros(NROIs, 2);
SubPixelOffsets = PixelOffsets;
for nn = 1:NROIs
    RowIndices = BlockROIs(nn, 1):BlockROIs(nn, 3);
    ColIndices = BlockROIs(nn, 2):BlockROIs(nn, 4);
    [Offset, SubOffset] = MIC_Reg3DTrans.findStackOffset(...
        Image1(RowIndices, ColIndices), Image2(RowIndices, ColIndices), ...
        [ImageSize, 1], [], [], 0, UseGPU);
    PixelOffsets(nn, :) = Offset(1:2);
    SubPixelOffsets(nn, :) = SubOffset(1:2);
end


end