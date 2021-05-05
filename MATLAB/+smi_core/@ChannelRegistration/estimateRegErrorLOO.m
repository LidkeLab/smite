function [SquaredError] = estimateRegErrorLOO(...
    TransformationType, TransformationParams, ...
    MovingCoordinates, FixedCoordinates)
%estimateRegErrorLOO estimates registration error by leave-one-out analysis
% This method will estimate the registration error as the squared
% difference between the fiducial coordinates and the transformed
% coordinates, with the transform used being computed for all other
% fiducial coordinates but the current point.  That is, this method
% estimates the registration error at each fiducial (control) point by a
% leave-one-out (LOO) analysis.
%
% INPUTS:
%   TransformationType: See smi_core.ChannelRegistration property of the
%                       same name.
%   TransformationParams: Cell array of additional parameters needed based
%                         on the TransformationType. For example, for
%                         TransformationType = 'lwm', this is a cell array
%                         with one element: NNeighborPoints. For
%                         TransformationType = 'polynomial', the one
%                         parameter is the polynomial degree (see usage
%                         below).
%   MovingCoordinates: Coordinates of points in the second fiducial.
%                      (Nx2 numeric array)
%   FixedCoordinates: Coordinates of points in the first fiducial.
%                     (Nx2 numeric array)
%
% OUTPUTS:
%   SquaredError: Squared error computed between FixedCoordinates and
%                 transformed MovingCoordinates, where the transform is
%                 found using all other points but the current point.
%                 (Pixels)(Nx1 numeric array)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Loop through each fiducial coordinate, compute a transform that excludes
% that point, and then determine the squared error after applying that
% transform to the fiducial coordinate.
NCoordinates = size(MovingCoordinates, 1);
IndexArray = 1:NCoordinates;
SquaredError = zeros(NCoordinates, 1);
for nn = IndexArray
    % Recompute the transform from all points but the current point.
    IndexArrayCurrent = IndexArray(IndexArray ~= nn);
    MovingCoordsCurrent = MovingCoordinates(IndexArrayCurrent, :);
    FixedCoordinatesCurrent = FixedCoordinates(IndexArrayCurrent, :);
    switch TransformationType
        case 'lwm'
            RegistrationTransform = fitgeotrans(...
                MovingCoordsCurrent, ...
                FixedCoordinatesCurrent, ...
                TransformationType, TransformationParams{1});
        case 'polynomial'
            RegistrationTransform = fitgeotrans(...
                MovingCoordsCurrent, ...
                FixedCoordinatesCurrent, ...
                TransformationType, TransformationParams{1});
        otherwise
            RegistrationTransform = fitgeotrans(...
                MovingCoordsCurrent, ...
                FixedCoordinatesCurrent, ...
                TransformationType);
    end
    
    % Apply the transform to the excluded coordinate and save the
    % squared error.
    SquaredError(nn) = ...
        smi_core.ChannelRegistration.estimateRegistrationError(...
        RegistrationTransform, ...
        MovingCoordinates(nn, :), FixedCoordinates(nn, :));
end


end