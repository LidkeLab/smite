function [MovingCoordinates, FixedCoordinates] = transformCoords(...
    RegistrationTransform, MovingCoordinates, FixedCoordinates)
%transformCoords transforms a set of coordinates with the given transform.
% This method transforms the input coordinates, with the set being
% transformed depending on the type of transform.  The purpose of this
% method is to eliminate the guesswork needed when applying the transform
% (e.g., 'lwm' doesn't have a forward transform, so we have to do
% "reverse" the effect of that transform). 
%
% INPUTS:
%   RegistrationTransform: RegistrationTransform is a MATLAB tform object
%                          describing the transform to be visualized.
%                          (tform object)
%   MovingCoordinates: Coordinates of points in the second fiducial.
%                      (Nx2 numeric array)
%   FixedCoordinates: Coordinates of points in the first fiducial.
%                     (Nx2 numeric array)
%
% OUTPUTS:
%   MovingCoordinates: Coordinates of points in the second fiducial.  These
%                      coordinates may or may not be transformed, depending
%                      on the type of transform used.
%                      (Nx2 numeric array)
%   FixedCoordinates: Coordinates of points in the first fiducial.  These
%                     coordinates may or may not be transformed, depending
%                     on the type of transform used.
%                     (Nx2 numeric array)

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Determine the transformation type and proceed as appropriate.
if ismember(class(RegistrationTransform), ...
        {'images.geotrans.LocalWeightedMeanTransformation2D', ...
        'images.geotrans.PolynomialTransformation2D', ...
        'images.geotrans.PiecewiseLinearTransformation2D'})
    % Since these transforms only have the transformPointsInverse() method,
    % we'll need to "reverse" the transform at each point (i.e., use the
    % inverse transform, see how far points moved, and then move the
    % initial points the same amount in the opposite direction).
    MovingCoordinatesTransformed = ...
        smi_core.ChannelRegistration.transformCoordsDirect(...
        RegistrationTransform, MovingCoordinates);
    MovingCoordinates = 2*MovingCoordinates - MovingCoordinatesTransformed;
else
    MovingCoordinates = ...
        smi_core.ChannelRegistration.transformCoordsDirect(...
        RegistrationTransform, MovingCoordinates);
end


end