function [ImageStack] = shiftImage(ImageStack, Shift, UseGPU)
%shiftImage shifts an image by the provided shift.
% This method uses the Fourier shift theorem to apply the provided 'Shift'
% to the 'ImageStack'.
% INPUTS:
%   ImageStack: Image(s) to be shifted. (YSizexXSizexNImages float)
%   Shift: Shift to be applied to ImageStack. (2x1 or 3x1 float)(YxXxZ)
%   UseGPU: Flag indicating GPU should be used. (Default = false)
%
% OUTPUTS:
%   ImageStack: Shifted image stack. (YSizexXSizexNImages float)

% Created by:
%   David J. Schodt (Lidke Lab 2021)


% Set defaults.
if (~exist('UseGPU', 'var') || isempty(UseGPU))
    UseGPU = false;
end

% Validate inputs.
ImSize = size(ImageStack, 1:3);
Shift = padarray(Shift, max(0, sum(ImSize>1)-numel(Shift)), 'post');

% Fourier transform the image stack.
if UseGPU
    ImageStack = gpuArray(ImageStack);
end
ImageStackFT = fftshift(fftn(ImageStack));

% Apply the shift.
[Kx, Ky, Kz] = meshgrid(1:ImSize(2), 1:ImSize(1), 1:ImSize(3));
ShiftedImFT = ImageStackFT .* exp(-2*pi*1i*Shift(1)*Ky/ImSize(1)) ...
    .* exp(-2*pi*1i*Shift(2)*Kx/ImSize(2)) ...
    .* exp(-2*pi*1i*Shift(3)*Kz/ImSize(3));

% Mask components beyond the Nyquist frequency.
NyquistEllipse = ((Kx - ImSize(2)/2) / (ImSize(2)/2)).^2 ...
    + ((Ky - ImSize(1)/2) / (ImSize(1)/2)).^2 ...
    + ((Kz - ImSize(3)/2) / (ImSize(3)/2)).^2;
ShiftedImFT(NyquistEllipse > 1) = 0;

% Invert the image back to the spatial domain.
ImageStack = abs(ifftn(fftshift(ShiftedImFT)));
if UseGPU
    ImageStack = gather(ImageStack);
end


end