function [PlotAxes] = visualizeRegistrationResults(PlotAxes, ...
    RegistrationTransform, ...
    MovingCoordinates, FixedCoordinates, ...
    MovingImage, FixedImage)
%visualizeRegistrationResults shows registration results on the fiducials.
% This method will demonstrate the registration performance on the
% fiducials by plotting appropriately transformed overlays of the raw data
% and the resulting coordinates.
%
% NOTE: For some transforms, we have to treat FixedCoordinates as the
%       moving pair (forward transforms don't exist for some transforms).  
%
% NOTE: Default values can be used (when available) by setting the
%       corresponding input to [].  For example, 
%       visualizeRegistrationResults([], TForm, Coords1, Coords2, ...
%           Fiducial1, Fiducial2) will use a default value for PlotAxes.
%
% INPUTS:
%   PlotAxes: Axes in which we will make the plots. (Default = gca())
%   RegistrationTransform: RegistrationTransform is a MATLAB tform object
%                          describing the transform to be visualized.
%                          (tform object)
%   MovingCoordinates: Coordinates of points in the second fiducial.
%                      (Nx2 numeric array)
%   FixedCoordinates: Coordinates of points in the first fiducial.
%                     (Nx2 numeric array)
%   MovingImage: Image of the second fiducial. (MxP numeric array)
%   FixedImage: Image of the first fiducial. (MxP numeric array)
%
% OUTPUTS:
%   PlotAxes: A MATLAB axes handle for the axes in which the 
%             visualizations are contained. 

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Set default parameters if needed.
if (~exist('PlotAxes', 'var') || isempty(PlotAxes) ...
        || ~isvalid(PlotAxes))
    PlotAxes = gca();
end

% Transform the input Fiducial2.
MovingImageTransformed = smi_core.ChannelRegistration.transformImages(...
    RegistrationTransform, MovingImage);

% Perform a full scale histogram stretch on the fiducials.
FixedImage = (FixedImage-min(FixedImage(:))) ...
    ./ max(max(FixedImage-min(FixedImage(:))));
MovingImageTransformed = ...
    (MovingImageTransformed-min(MovingImageTransformed(:))) ...
    ./ max(max(MovingImageTransformed-min(MovingImageTransformed(:))));

% Create an RGB image of Fiducial1 and Fiducial2Transformed overlain.
FiducialOverlay = zeros([size(FixedImage), 3]);
FiducialOverlay(:, :, 2) = FixedImage;
FiducialOverlay(:, :, 1) = MovingImageTransformed;
FiducialOverlay(:, :, 3) = MovingImageTransformed;

% Transform the coordinates.
[MovingCoordinatesTrans, FixedCoordinatesTrans] = ...
    smi_core.ChannelRegistration.transformCoords(...
    RegistrationTransform, MovingCoordinates, FixedCoordinates);

% Plot the fiducial images and the coordinates.
ImshowCorrection = [0.5, 0.5];
imshow(FiducialOverlay, [], 'Parent', PlotAxes, ...
    'XData', [1, size(FiducialOverlay, 2)]+ImshowCorrection, ...
    'YData', [1, size(FiducialOverlay, 1)]+ImshowCorrection)
line(PlotAxes, ...
    FixedCoordinatesTrans(:, 1), FixedCoordinatesTrans(:, 2), ...
    'Marker', 'o', 'Color', [0, 1, 0], ...
    'LineStyle', 'none', 'LineWidth', 1.5)
line(PlotAxes, ...
    MovingCoordinatesTrans(:, 1), MovingCoordinatesTrans(:, 2), ...
    'Marker', '+', 'Color', [1, 0, 1], ...
    'LineStyle', 'none', 'LineWidth', 1.5)
title(PlotAxes, ...
    {'Fixed fiducial - green', 'Moving fiducial - magenta'})


end