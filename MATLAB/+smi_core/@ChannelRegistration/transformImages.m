function [TransformedImages] = transformImages(...
    RegistrationTransform, Images)
%transformImages transforms a set of images with the given transform.
% This method will transform the images given in Images using the
% transform in RegistrationTransform.
%
% NOTE: This method is just a wrapper for the MATLAB methods imref2D() and
%       imwarp().
%
% INPUTS:
%   RegistrationTransform: A MATLAB tform object containing information
%                          about the transformation to be used
%                          (tform object)
%   Images: A stack of images to be transformed (MxNxP numeric array)
%
% OUTPUTS:
%   TransformedImages: Input Images transformed using RegistrationTransform
%                      (MxNxP numeric array)

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Apply the transform to the images.
if ~isempty(RegistrationTransform)
    ImRef2DObject = imref2d(size(Images(:, :, 1)));
    TransformedImages = imwarp(Images, RegistrationTransform, ...
        'OutputView', ImRef2DObject);
end


end