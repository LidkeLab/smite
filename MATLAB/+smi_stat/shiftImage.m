function [ImageStack, Params] = shiftImage(ImageStack, Shift, Params)
%shiftImage shifts an image by the provided shift.
% This method uses the Fourier shift theorem to apply the provided 'Shift'
% to the 'ImageStack'.
% INPUTS:
%   ImageStack: Image(s) to be shifted. (YSizexXSizexNImages float)
%   Shift: Shift to be applied to ImageStack. (2x1 or 3x1 float)(YxXxZ)
%   Params: Structure of parameters.
%           UseGPU: Flag indicating GPU should be used. (Default = false)
%           RemoveEdges: Flag indicating false edges should be removed.
%                        For example, if we shift an image by 
%                        [0.3; -2.2; 1.6],
%                        we'll set row ceil(0.3)=1 to PadValue, column
%                        X-ceil(2.2)=X-3 to PadValue, and images
%                        1:ceil(1.6)=1:2 to PadValue. (Default = false)
%           PadValue: Value used to replace the removed periodic edges when
%                     RemoveEdges = true.  If set to [], PadValue will be
%                     set to mean(ImageStack(:)). (Default = [])
%
% OUTPUTS:
%   ImageStack: Shifted image stack. (YSizexXSizexNImages float)
%   Params: Input parameters padded with defaults.

% Created by:
%   David J. Schodt (Lidke Lab 2021)


% Set defaults.
DefaultParams.UseGPU = false;
DefaultParams.RemoveEdges = false;
DefaultParams.PadValue = [];
if (~exist('Params', 'var') || isempty(Params))
    Params = DefaultParams;
else
    Params = smi_helpers.padStruct(Params, DefaultParams);
end

% Validate inputs.
ImSize = size(ImageStack, 1:3);
Shift = padarray(Shift, max(0, sum(ImSize>1)-numel(Shift)), 'post');

% Fourier transform the image stack.
if Params.UseGPU
    ImageStack = gpuArray(ImageStack);
end
ImageStackFT = fftshift(fftn(ImageStack));

% Apply the shift.
[Kx, Ky, Kz] = meshgrid(1:ImSize(2), 1:ImSize(1), 1:ImSize(3));
ShiftedImFT = ImageStackFT .* exp(-2*pi*1i*Shift(1)*Ky/ImSize(1)) ...
    .* exp(-2*pi*1i*Shift(2)*Kx/ImSize(2)) ...
    .* exp(-2*pi*1i*Shift(3)*Kz/ImSize(3));

% Mask components beyond the Nyquist frequency.
NyquistEllipse = ((Kx - ImSize(2)/2 - 0.5) / (ImSize(2)/2 - 0.5)).^2 ...
    + ((Ky - ImSize(1)/2 - 0.5) / (ImSize(1)/2 - 0.5)).^2 ...
    + ((Kz - ImSize(3)/2 - 0.5) / (ImSize(3)/2 - 0.5)).^2;
ShiftedImFT(NyquistEllipse > 1) = 0;

% Invert the image back to the spatial domain.
ImageStack = abs(ifftn(fftshift(ShiftedImFT)));

% Remove the periodic edges if requested.
if Params.RemoveEdges
    if isempty(Params.PadValue)
        Params.PadValue = mean(ImageStack(:));
    end
    YDelete = (1:min(ImSize(1), ceil(abs(Shift(1))))) ...
        + (Shift(1)<0)*mod(floor(Shift(1)), ImSize(1));
    ImageStack(YDelete, :, :) = Params.PadValue;
    XDelete = (1:min(ImSize(2), ceil(abs(Shift(2))))) ...
        + (Shift(2)<0)*mod(floor(Shift(2)), ImSize(2));
    ImageStack(:, XDelete, :) = Params.PadValue;
    ZDelete = (1:min(ImSize(3), ceil(abs(Shift(3))))) ...
        + (Shift(3)<0)*mod(floor(Shift(3)), ImSize(3));
    ImageStack(:, :, ZDelete) = Params.PadValue;
end
if Params.UseGPU
    ImageStack = gather(ImageStack);
end


end