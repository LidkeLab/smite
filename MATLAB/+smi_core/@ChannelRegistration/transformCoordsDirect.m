function [TransformedCoordinates] = transformCoordsDirect(...
    RegistrationTransform, Coordinates)
%transformCoordsDirect transforms a set of coordinates.
% This method transforms the input coordinates directly.  For some
% transforms, the inverse transform is applied (nothing in the output will
% distinguish this).  For this reason, it is recommended that you use
% transformCoords() instead of this method.
%
% INPUTS:
%   RegistrationTransform: RegistrationTransform is a MATLAB tform object
%                          describing the transform to be visualized.
%                          (tform object)
%   Coordinates: Coordinates to be transformed. (Nx2 numeric array)
%
% OUTPUTS:
%   TransformedCoordinates: Input Coordinates transformed by
%                           RegistrationTransform. (Nx2 numeric array)

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Determine the transformation type and proceed as appropriate.
if ismember(class(RegistrationTransform), ...
        {'images.geotrans.LocalWeightedMeanTransformation2D', ...
        'images.geotrans.PolynomialTransformation2D', ...
        'images.geotrans.PiecewiseLinearTransformation2D'})
    % None of these transforms have a transformPointsForward()
    % implementation, so we must use the inverse on our reference points
    % (We could just compute the transform treating the reference as
    % moving, but then the transform can't be applied to our raw data since
    % MATLAB doesn't have an inverse of imwarp().  For SPT, it's nice to
    % have transformed raw data for use in movies.)
    [TransformedCoordinates(:, 1), TransformedCoordinates(:, 2)] = ...
        transformPointsInverse(RegistrationTransform, ...
        Coordinates(:, 1), Coordinates(:, 2));
else
    [TransformedCoordinates(:, 1), TransformedCoordinates(:, 2)] = ...
        transformPointsForward(RegistrationTransform, ...
        Coordinates(:, 1), Coordinates(:, 2));
end


end