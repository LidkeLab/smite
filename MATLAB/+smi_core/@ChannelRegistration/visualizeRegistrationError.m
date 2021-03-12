function [PlotAxes] = visualizeRegistrationError(PlotAxes, ...
    RegistrationTransform, MovingCoordinates, FixedCoordinates, ...
    FOV, GridSpacing)
%visualizeRegistrationError visualizes the error in RegistrationTransform
% This method will plot the error of RegistrationTransform.  This is done
% by finding the registration error as the root mean square difference
% between the fiducial coordinates being compared after transforming (e.g.,
% fiducial 1 localizations and transformed fiducial 2 localizations) and
% interpolating between the points to produce a smooth heat map type image
% estimate of the error. 
%
% NOTE: For some transforms, we have to treat FixedCoordinates as the
%       moving pair (forward transforms don't exist for some transforms).  
%
% NOTE: Default values can be used (when available) by setting the
%       corresponding input to [].  For example, 
%       visualizeRegistrationError([], TForm, Coords1, Coords2) will use a
%       default value for PlotAxes.
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
%   FOV: Field of view in which we will estimate the error.
%        (Pixels)(1x4, [XStart, YStart, XEnd, YEnd])
%        (Default is set to encompass all points in FixedCoordinates)
%   GridSpacing: Spacing of a grid made across FrameSize.  The grid will
%                consist of points
%                [0:GridSpacing:FrameSize(1), 0:GridSpacing:FrameSize(2)] 
%                at which we will apply RegistrationTransform.
%                (Pixels)(Default = 1)
% OUTPUTS:
%   PlotAxes: A MATLAB axes handle for the axes in which the 
%             visualizations are contained. 

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Set default parameters and modify inputs.
FixedCoordinates = double(FixedCoordinates);
MovingCoordinates = double(MovingCoordinates);
if (~exist('PlotAxes', 'var') || isempty(PlotAxes) ...
        || ~isvalid(PlotAxes))
    PlotAxes = gca();
end
if (~exist('GridSpacing', 'var') || isempty(GridSpacing))
    GridSpacing = 1;
end
if (~exist('FOV', 'var') || isempty(FOV))
    FOV = [floor(min(FixedCoordinates(:, 1))-GridSpacing), ...
        floor(min(FixedCoordinates(:, 2))-GridSpacing), ...
        ceil(max(FixedCoordinates(:, 1))+GridSpacing), ...
        ceil(max(FixedCoordinates(:, 2))+GridSpacing)];
end

% Compute the squared difference between FixedCoordinates and transformed 
% MovingCoordinates.
SquaredError = smi_core.ChannelRegistration.estimateRegistrationError(...
    RegistrationTransform, MovingCoordinates, FixedCoordinates);

% Create an interpolant function handle for our registration error.
InterpolantFunction = scatteredInterpolant(...
    FixedCoordinates(:, 1), FixedCoordinates(:, 2), ...
    double(SquaredError), ...
    'natural', 'linear');

% Estimate the squared error on a grid.
[XGrid, YGrid] = ...
    meshgrid(FOV(1):GridSpacing:FOV(3), FOV(2):GridSpacing:FOV(4));
SquaredErrorGrid = InterpolantFunction(XGrid(:), YGrid(:));

% Display the grid of estimated squared errors.
SquaredErrorGrid = reshape(SquaredErrorGrid, size(XGrid));
ImshowCorrection = [0.5, -0.5];
imshow(SquaredErrorGrid, [], 'Parent', PlotAxes, ...
    'XData', FOV([1, 3])+ImshowCorrection, ...
    'YData', FOV([2, 4])+ImshowCorrection)
colormap(PlotAxes, 'parula')
title(PlotAxes, ...
    'Squared registration error (interpolated) in full field of view')
xlabel(PlotAxes, 'X (pixels)')
ylabel(PlotAxes, 'Y (pixels)')
axis(PlotAxes, 'tight')
ColorBar = colorbar(PlotAxes);
ColorBar.Label.String = 'Squared error (pixels^2)';

% Plot the true points as a visual reference.
hold(PlotAxes, 'on');
line(PlotAxes, FixedCoordinates(:, 1), FixedCoordinates(:, 2), ...
    'Color', 'k', 'Marker', '.', 'LineStyle', 'none')


end