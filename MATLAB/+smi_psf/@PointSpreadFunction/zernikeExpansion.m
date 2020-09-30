function [NollCoef,ImageZ]=zernikeExpansion(Image,ZStruct)
%zernikeExpansion Expand Image into Zernike Moments
%
%
%   Noll ordering is used as a linear index.
%
%   More about formula and conventions.
%
% INPUTS:
%   Image:          Image
%   ZStruct:        Pre-computed images for quick expansion and sum
%
% OUTPUTS:
%   NollCoefs:      Noll Coefficients
%   ImageZ:         Zernike Polynomial Image
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%

Image=Image(:);

NollCoef=(ZStruct.ZImages*Image).*ZStruct.Norms';   

if nargout>1 %Create Image
    [ImageZ]=smi_psf.PointSpreadFunction.zernikeSum(NollCoef,ZStruct);    
end
