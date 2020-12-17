function [SquaredError] = estimateRegistrationError(...
    RegistrationTransform, Coords1, Coords2)
%estimateRegistrationError estimates the registration error.
% This method will estimate the registration error as the squared
% difference between the fiducial coordinates being compared after 
% transforming (e.g., fiducial 1 localizations and transformed fiducial 2
% localizations).
%
% NOTE: If you need a single number to quantify your registration error,
%       a good estimate would be RMSE = sqrt(mean(SquaredError))
%
% NOTE: Coords2 will always be the set of coordinates that will be
%       transformed.
%
% INPUTS:
%   RegistrationTransform: RegistrationTransform is a MATLAB tform object
%                          describing the transform to be visualized.
%                          (tform object)
%   Coords1: Coordinates of points in the first fiducial.
%            (Pixels)(Nx2 numeric array)
%   Coords2: Coordinates of points in the second fiducial that will be
%            transformed by RegistrationTransform. 
%            (Pixels)(Nx2 numeric array)
%
% OUTPUTS:
%   SquaredError: Squared error computed between Coords1 and transformed
%                 Coords2. This is returned for convenience in other 
%                 codes (i.e., instead of the RMSE)
%                 (Pixels)(Nx1 numeric array)

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Transform the input Coords2 and compute the error.
Coords2Transformed = smi_core.ChannelRegistration.transformCoords(...
    RegistrationTransform, Coords2);
SquaredError = sum((Coords1 - Coords2Transformed).^2, 2);


end