function [Image]=zernikeSum(NollCoefs,ZStruct)
%zernikeSum Generate an image from Noll Coefficients,
%
%
%   Noll ordering is used as a linear index.
%
%   More about formula and conventions.
%
% INPUTS:
%   NollCoefs:      Noll Coefficients (Nx1);
%   ZStruct         Pre-computed images for quick expansion and sum
%
% OUTPUTS:
%   Image:          Zernike Polynomial Image
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%

if length(NollCoefs)>ZStruct.NMax
    error('zernikeSum:: ZStruct too small')
end

if length(NollCoefs)<ZStruct.NMax
    ZImages=ZStruct.ZImages(1:length(NollCoefs),:);
else
    ZImages=ZStruct.ZImages;
end

NollCoefs=diag(NollCoefs);
Image=sum(NollCoefs*ZImages,1);
Image=reshape(Image,[ZStruct.SZ,ZStruct.SZ]);






