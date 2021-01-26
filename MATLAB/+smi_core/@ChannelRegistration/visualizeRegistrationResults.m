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
%           Fiducial1, Fiducial2) will use a default value for PlotFigure.
%
% INPUTS:
%   PlotFigure: Figure in which we will make the plots. (Default = gcf())
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
%   PlotFigure: A MATLAB figure handle for the figure in which the 
%               visualizations are contained. 

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Set default parameters if needed.
if (~exist('PlotFigure', 'var') || isempty(PlotFigure) ...
        || ~isvalid(PlotFigure))
    PlotFigure = gcf();
end

% Transform the moving fiducial..
MovingImageTransformed = smi_core.ChannelRegistration.transformImages(...
    RegistrationTransform, MovingImage);

% Perform a full scale histogram stretch on the fiducials.
FixedImage = (FixedImage-min(FixedImage(:))) ...
    ./ max(max(FixedImage-min(FixedImage(:))));
MovingImage = (MovingImage-min(MovingImage(:))) ...
    ./ max(max(MovingImage-min(MovingImage(:))));
MovingImageTransformed = ...
    (MovingImageTransformed-min(MovingImageTransformed(:))) ...
    ./ max(max(MovingImageTransformed-min(MovingImageTransformed(:))));

% Create an RGB image of the fixed and moving fiducials (both before and 
% after the transform).
FiducialOverlay = zeros([size(FixedImage), 3]);
FiducialOverlay(:, :, 2) = FixedImage;
FiducialOverlay(:, :, 1) = MovingImage;
FiducialOverlay(:, :, 3) = MovingImage;
FiducialOverlayTrans = zeros([size(FixedImage), 3]);
FiducialOverlayTrans(:, :, 2) = FixedImage;
FiducialOverlayTrans(:, :, 1) = MovingImageTransformed;
FiducialOverlayTrans(:, :, 3) = MovingImageTransformed;

% Transform the coordinates.
[MovingCoordinatesTrans, FixedCoordinatesTrans] = ...
    smi_core.ChannelRegistration.transformCoords(...
    RegistrationTransform, MovingCoordinates, FixedCoordinates);

% Compute the RMSE.
RMSEBefore = sqrt(mean(sum((MovingCoordinates - FixedCoordinates).^2, 2)));
RMSEAfter = sqrt(mean(...
    smi_core.ChannelRegistration.estimateRegistrationError(...
    RegistrationTransform, MovingCoordinates, FixedCoordinates)));

% Plot the fiducial images and the coordinates.
ImshowCorrection = [0.5, 0.5];
PlotAxesLeft = subplot(1, 2, 1, 'Parent', PlotFigure);
imshow(FiducialOverlay, [], 'Parent', PlotAxesLeft, ...
    'XData', [1, size(FiducialOverlayTrans, 2)]+ImshowCorrection, ...
    'YData', [1, size(FiducialOverlayTrans, 1)]+ImshowCorrection)
line(PlotAxesLeft, ...
    FixedCoordinates(:, 1), FixedCoordinates(:, 2), ...
    'Marker', 'x', 'MarkerSize', 8, 'Color', [0, 1, 0], ...
    'LineStyle', 'none', 'LineWidth', 1.5)
line(PlotAxesLeft, ...
    MovingCoordinates(:, 1), MovingCoordinates(:, 2), ...
    'Marker', '+', 'MarkerSize', 8, 'Color', [1, 0, 1], ...
    'LineStyle', 'none', 'LineWidth', 1.5)
title(PlotAxesLeft, ...
    {sprintf('BEFORE: RMSE = %.3f pixels', RMSEBefore), ...
    'Fixed fiducial - green', 'Moving fiducial - magenta'})
xlabel(PlotAxesLeft, 'X (pixels)')
ylabel(PlotAxesLeft, 'Y (pixels)')
PlotAxesRight = subplot(1, 2, 2, 'Parent', PlotFigure);
imshow(FiducialOverlayTrans, [], 'Parent', PlotAxesRight, ...
    'XData', [1, size(FiducialOverlayTrans, 2)]+ImshowCorrection, ...
    'YData', [1, size(FiducialOverlayTrans, 1)]+ImshowCorrection)
line(PlotAxesRight, ...
    FixedCoordinatesTrans(:, 1), FixedCoordinatesTrans(:, 2), ...
    'Marker', 'x', 'MarkerSize', 8, 'Color', [0, 1, 0], ...
    'LineStyle', 'none', 'LineWidth', 1.5)
line(PlotAxesRight, ...
    MovingCoordinatesTrans(:, 1), MovingCoordinatesTrans(:, 2), ...
    'Marker', '+', 'MarkerSize', 8, 'Color', [1, 0, 1], ...
    'LineStyle', 'none', 'LineWidth', 1.5)
title(PlotAxesRight, ...
    {sprintf('AFTER: RMSE = %.3f pixels', RMSEAfter), ...
    'Fixed fiducial - green', 'Moving fiducial - magenta'})
xlabel(PlotAxesRight, 'X (pixels)')
ylabel(PlotAxesRight, 'Y (pixels)')

end