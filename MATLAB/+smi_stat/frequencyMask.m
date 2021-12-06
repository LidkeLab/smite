function [FreqMask, FreqSqEllipse, YMesh, XMesh, ZMesh] = ...
    frequencyMask(ImSize, FreqCutoff)
%frequencyMask prepares a boolean mask defining a frequency cutoff.
% This method prepares a boolean array that defines the indices of a
% Fourier transform image that are <= or > the cutoff frequency.  For
% Fourier transforms with radially increasing frequency (i.e., max
% frequency at center of image), the output mask defines signal below or at
% the cutoff frequency.  For Fourier transforms with radially decreasing
% frequency, the mask defines signal above the cutoff frequency.
%
% INPUTS:
%   ImSize: Size of the Fourier transform image, obtained by, e.g.,
%           size(FT). (RowsxColsxZ)
%   FreqCutoff: Cutoff (cuton) frequency of the mask (see note above).
%               (1/pixels)(Default = Nyquist frequency)           
%
% OUTPUTS:
%   FreqMask: Boolean array defining the indices to be masked. 
%             (1xprod(ImSize))
%   FreqSqEllipse: Ellipse defining the squared image frequencies (can be 
%                  restored to as an image by doing reshape(FreqEllipse, 
%                  ImSize)). This is defined by the convention that 
%                  frequency is radially increasing from the image center, 
%                  with each pixel defining the frequency at the pixel 
%                  center (as opposed to, e.g., the far edge of the pixel, 
%                  which would give a more restrictive mask). 
%                  (1xprod(ImSize))
%   YMesh, XMesh, ZMesh: Results from 
%                        meshgrid(1:ImSize(2), 1:ImSize(1), 1:ImSize(3)).
%                        This is the most expensive piece of this code, so
%                        it's returned just in case the user still needs
%                        it.

% Created by:
%   David J. Schodt (Lidke Lab 2021)


% Prepare defaults.
FNyquist = 0.5;
if (~exist('FreqCutoff', 'var') || isempty(FreqCutoff))
    % Define the default frequency as the Nyquist frequency 
    % FNyquist = FSampling / 2 = (1 sample/pixel)/2 = 0.5 / pixel
    FreqCutoff = FNyquist;
end

% Reshape and pad inputs if needed.
if isrow(ImSize)
    ImSize = ImSize.';
end
ImSize = padarray(ImSize, [3-numel(ImSize), 0], 1, 'post');

% Compute the image frequencies at each pixel of an image sized 'ImSize'.
[XMesh, YMesh, ZMesh] = meshgrid(1:ImSize(2), 1:ImSize(1), 1:ImSize(3));
FreqSqEllipse = ((YMesh-ImSize(1)/2-0.5) / (ImSize(1)/2)).^2 ...
    + ((XMesh-ImSize(2)/2-0.5) / (ImSize(2)/2)).^2 ...
    + ((ZMesh-ImSize(3)/2-0.5) / (ImSize(3)/2)).^2;

% Compute the binary mask.
% NOTE: FreqEllipse==1 lies right on the Nyquist frequency of the image,
%       so we can compare the ellipse to FreqCutoff * (RNyquist/FNyquist) =
%       FreqCutoff * (1/FNyquist) to get the desired cutoff.
FreqMask = (FreqSqEllipse > (FreqCutoff*(1/FNyquist))^2);


end