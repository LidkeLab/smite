function [DividedImages, ImageROIs] = subdivideImage(Image, SubROISize)
%subdivideImage divvys up an image into sub-ROIs of the image.
% This method divvys up the input image Image into sub-ROIs matching
% the size specified in SubROISize, with edge cases taken to be the largest
% ROI possible (that is smaller than the nominal sub-ROI size).
%
% INPUT:
%   Image: The image which will be roughly divided into SubROISize sized
%          images. (MxN numeric array)
%   SubROISize: The nominal size of the output subdivided images.
%               (Pixels)(2x1 array)(Default = size(Image))
%
% OUTPUT:
%   DividedImages: Cell array of the images subdivided from Image.
%                  (NROIsx1 cell array)
%   ImageROIs: ROIs of the regions corresponding to the divided images.
%              (NROIsx4 array)([YStart, XStart, YEnd, XEnd])
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
ImageSize = size(Image, [1, 2]);
if (~exist('SubROISize', 'var') || isempty(SubROISize))
    SubROISize = ImageSize.';
end

% Ensure arrays are properly shaped.
if isrow(SubROISize)
    SubROISize = SubROISize.';
end

% Define the ROIs of the subdivided images.
NDivisions = ceil(ImageSize.' ./ SubROISize);
YStart = repmat(1 + SubROISize(1)*(0:(NDivisions(1)-1)).', ...
    [NDivisions(2), 1]);
XStart = repelem(1 + SubROISize(2)*(0:(NDivisions(2)-1)).', ...
    NDivisions(1));
ImageROIs = [YStart, XStart, ...
    min(ImageSize(1), YStart+SubROISize(1)-1), ...
    min(ImageSize(2), XStart+SubROISize(2)-1)];

% Loop through each ROI and isolate that section of the image.
NROIs = prod(NDivisions);
DividedImages = cell(NROIs, 1);
for nn = 1:NROIs
    DividedImages{nn} = Image(ImageROIs(nn, 1):ImageROIs(nn, 3), ...
        ImageROIs(nn, 2):ImageROIs(nn, 4), :);
end


end