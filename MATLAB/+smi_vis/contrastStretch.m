function [Image] = contrastStretch(Image, MinMax, ...
    PercentileCeiling, MinScaleIntensity)
%contrastStretch scales images to the range MinMax.
% This method performs a full scale histogram stretch of 'Image' such that
% the stretched pixel values lie in the range [MinMax(1), MinMax(2)].
%
% INPUTS:
%   Image: Array of pixel values that will be stretched. (float array)
%   MinMax: Array containing the minimum and maximum pixel value after
%           scaling. (Default = [0, 1])
%   MinScaleIntensity: Minimum scaling intensity (useful for noisy data,
%                      so that the scaling doesn't just brighten the noise)
%                      (Default = 1)
%   PercentileCeiling: Percentile ceiling of pixel values in the raw data
%                      above which values are clipped. (Default = 100)
%
% OUTPUTS:
%   Image: Input array 'Image' scaled to the range 'MinMax'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults and check some inputs.
if (~exist('MinMax', 'var') || isempty(MinMax))
    MinMax = [0; 1];
end
MinMax = sort(MinMax);
if (~exist('MinScaleIntensity', 'var') || isempty(MinScaleIntensity))
    MinScaleIntensity = 1;
end
if (~exist('PercentileCeiling', 'var') || isempty(PercentileCeiling))
    PercentileCeiling = 100;
end

% If the image is identically zero, return.
if ~any(Image(:))
    Image = Image + MinMax(1);
    return
end

% Scale the 'Image' array.
if (PercentileCeiling ~= 100)
    MaxIntensity = ...
        max(prctile(Image(:), PercentileCeiling), MinScaleIntensity);
    Image(Image > MaxIntensity) = MaxIntensity;
end
Image = (Image-min(Image(:))) ...
    * (MinMax(2)-MinMax(1))/(max(Image(:))-min(Image(:))) ...
    + MinMax(1);


end