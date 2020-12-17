function [TransformedImages] = transformImages(Images, ...
    RegistrationTransform)
%transformImages transforms a set of images with the given transform.
% This method will transform the images given in Images using the
% transform in RegistrationTransform.
%
% NOTE: This method is just a wrapper for the MATLAB methods imref2D() and
%       imwarp().
%
% INPUTS:
%   Images: A stack of images to be transformed (MxNxP numeric array)
%   RegistrationTransform: A MATLAB tform object containing information
%                          about the transformation to be used 
%                          (tform object)
%
% OUTPUTS: 
%   TransformedImages: Input Images transformed usingRegistrationTransform
%                      (MxNxP numeric array)

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Apply the transform to the images.
ImRef2DObject = imref2d(size(Images(:, :, 1)));
TransformedImages = imwarp(Images, RegistrationTransform, ...
    'OutputView', ImRef2DObject);


end