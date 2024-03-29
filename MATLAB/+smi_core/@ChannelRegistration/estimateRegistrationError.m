function [SquaredError] = estimateRegistrationError(...
    RegistrationTransform, MovingCoordinates, FixedCoordinates)
%estimateRegistrationError estimates the registration error.
% This method will estimate the registration error as the squared
% difference between the fiducial coordinates being compared after 
% transforming (e.g., fiducial 1 localizations and transformed fiducial 2
% localizations).
%
% NOTE: If you need a single number to quantify your registration error,
%       a good estimate would be RMSE = sqrt(mean(SquaredError))
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
%   SquaredError: Squared error computed between FixedCoordinates and
%                 transformed MovingCoordinates. This is returned for 
%                 convenience in other codes (i.e., instead of the RMSE)
%                 (Pixels)(Nx1 numeric array)

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Transform the input Coords2 and compute the error.
if isempty(RegistrationTransform)
    SquaredError = NaN(size(MovingCoordinates, 1), 1);
    return
end
MovingCoordinates = ...
    smi_core.ChannelRegistration.transformCoords(...
    RegistrationTransform, MovingCoordinates);
SquaredError = sum((FixedCoordinates - MovingCoordinates).^2, 2);


end