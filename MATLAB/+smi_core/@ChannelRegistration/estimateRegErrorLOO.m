function [SquaredError] = estimateRegErrorLOO(obj)
%estimateRegErrorLOO estimates registration error by leave-one-out analysis
% This method will estimate the registration error as the squared
% difference between the fiducial coordinates and the transformed
% coordinates, with the transform used being computed for all other
% fiducial coordinates but the current point.  That is, this method
% estimates the registration error at each fiducial (control) point by a
% leave-one-out (LOO) analysis.  This method should typically be called
% after/at the end of running obj.findTransform(), as several pieces of
% code in that method will populate class properties needed in this
% analysis.
%
% OUTPUTS:
%   SquaredError: Squared error computed between the fiducial coordinates
%                 and the transformed fiducial coordinates.  Each entry of
%                 the cell array corresponds with the LOO squared error at
%                 each coordinate for each of the NFiducials.  That is,
%                 SquaredError{ff} has the LOO squared error at each
%                 coordinate for fiducial ff transforming to fiducial 1.
%                 (Pixels)(NFiducialsx1 cell array)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Loop through each fiducial coordinate, compute a transform that excludes
% that point, and then determine the squared error after applying that
% transform to the fiducial coordinate.
NFiducials = numel(obj.Coordinates);
SquaredError = cell(NFiducials, 1);
SquaredError{1} = zeros(size(obj.Coordinates{1}, 1), 1);
for ff = 2:NFiducials
    NCoordinates = size(obj.Coordinates{ff}, 1);
    IndexArray = 1:NCoordinates;
    SquaredError{ff} = zeros(NCoordinates, 1);
    for nn = IndexArray
        % Recompute the transform from all points but the current point.
        IndexArrayCurrent = IndexArray(IndexArray ~= nn);
        MovingCoordsCurrent = ...
            obj.Coordinates{ff}(IndexArrayCurrent, :, 2);
        ReferenceCoordsCurrent = ...
            obj.Coordinates{ff}(IndexArrayCurrent, :, 1);
        switch obj.TransformationType
            case 'lwm'
                RegistrationTransform = fitgeotrans(...
                    MovingCoordsCurrent, ...
                    ReferenceCoordsCurrent, ...
                    obj.TransformationType, obj.NNeighborPoints);
            case 'polynomial'
                RegistrationTransform = fitgeotrans(...
                    MovingCoordsCurrent, ...
                    ReferenceCoordsCurrent, ...
                    obj.TransformationType, obj.PolynomialDegree);
            otherwise
                RegistrationTransform = fitgeotrans(...
                    MovingCoordsCurrent, ...
                    ReferenceCoordsCurrent, ...
                    obj.TransformationType);
        end
        
        % Apply the transform the the excluded coordinate and save the
        % squared error.
        SquaredError{ff}(nn, 1) = sum((obj.Coordinates{ff}(nn, :, 1) ...
            - obj.transformCoords(RegistrationTransform, ...
            obj.Coordinates{ff}(nn, :, 2))).^2);
    end
end


end