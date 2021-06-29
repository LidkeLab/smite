function [Image] = contrastStretch(Image, Low, High)
%contrastStretch scales images to the range [Low, High].
% This method performs a full scale histogram stretch of 'Image' such that
% the stretched pixel values lie in the range [Low, High].
%
% INPUTS:
%   Image: Array of pixel values that will be stretched. (float array)
%   Low: Minimum pixel value in the output array 'Image'.
%        (numeric scalar)(Default = 0)
%   High: Maximum pixel value in the output array 'Image'.
%         (numeric scalar)(Default = 1)
%
% OUTPUTS:
%   Image: Input array 'Image' scaled to the range [Low, High].

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults and check some inputs.
if (~exist('Low', 'var') || isempty(Low))
    Low = 0;
end
if (~exist('High', 'var') || isempty(High))
    High = 1;
end
assert(Low < High, ...
    'contrastStretch(): Input ''Low'' must be less than ''High''')

% Scale the 'Image' array.
Image = (Image-min(Image(:))) * (High-Low)/(max(Image(:))-min(Image(:))) ...
    + Low;


end