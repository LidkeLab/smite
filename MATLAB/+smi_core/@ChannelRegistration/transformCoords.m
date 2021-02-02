function [MovingCoordsTransformed] = ...
    transformCoords(RegistrationTransform, MovingCoords)
%transformCoords transforms a set of coordinates with the given transform.
% This method transforms the input coordinates, with the set being
% transformed depending on the type of transform.  The purpose of this
% method is to overcome the annoyances/issues that may arise when applying 
% the transform directly (e.g., 'lwm' doesn't have a forward transform, so
% we have to do "reverse" the effect of that transform). 
%
% INPUTS:
%   RegistrationTransform: RegistrationTransform is a MATLAB tform object
%                          describing the transform to be visualized.
%                          (tform object)
%   MovingCoords: Coordinates of points in the "moving" fiducial.
%                 (Nx2 numeric array)
%
% OUTPUTS:
%   MovingCoordsTransformed: Input 'MovingCoords' transformed by 
%                            'RegistrationTransform'. 
%                            (Nx2 numeric array)

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
    MovingCoordsTransformed = ...
        smi_core.ChannelRegistration.transformCoordsDirect(...
        RegistrationTransform, MovingCoords);
    MovingCoordsTransformed = 2*MovingCoords ...
        - MovingCoordsTransformed;
else
    MovingCoordsTransformed = ...
        smi_core.ChannelRegistration.transformCoordsDirect(...
        RegistrationTransform, MovingCoords);
end


end