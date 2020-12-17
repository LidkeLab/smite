function [PlotAxes] = visualizeRegistrationError(PlotAxes, ...
    RegistrationTransform, Coords1, Coords2, FOV, GridSpacing)
%visualizeRegistrationError visualizes the error in RegistrationTransform
% This method will plot the error of RegistrationTransform.  This is done
% by finding the registration error as the root mean square difference
% between the fiducial coordinates being compared after transforming (e.g.,
% fiducial 1 localizations and transformed fiducial 2 localizations) and
% interpolating between the points to produce a smooth heat map type image
% estimate of the error. 
%
% NOTE: Coords2 will always be the set of coordinates that will be
%       transformed.
%
% NOTE: Default values can be used (when available) by setting the
%       corresponding input to [].  For example, 
%       visualizeCoordTransform([], TForm, FrameSize) will use a default
%       value for PlotAxes.
%
% INPUTS:
%   PlotAxes: Axes in which we will make the plots. (Default = gca())
%   RegistrationTransform: RegistrationTransform is a MATLAB tform object
%                          describing the transform to be visualized.
%                          (tform object)
%   Coords1: Coordinates of points in the first fiducial.
%            (Nx2 numeric array)
%   Coords2: Coordinates of points in the second fiducial that will be
%            transformed by RegistrationTransform. (Nx2 numeric array)
%   FOV: Field of view in which we will estimate the error.
%        (Pixels)(1x4, [XStart, YStart, XEnd, YEnd])
%        (Default set to encompass all points in Coords1)
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
Coords1 = double(Coords1);
Coords2 = double(Coords2);
if (~exist('PlotAxes', 'var') || isempty(PlotAxes) ...
        || ~isvalid(PlotAxes))
    PlotAxes = gca();
end
if (~exist('FOV', 'var') || isempty(FOV))
    FOV = [floor(min(Coords1(:, 2))), floor(min(Coords1(:, 1))), ...
        ceil(max(Coords1(:, 2))), ceil(max(Coords1(:, 1)))];
end
if (~exist('GridSpacing', 'var') || isempty(GridSpacing))
    GridSpacing = 1;
end

% Compute the squared difference between Coords1 and transformed Coords2.
SquaredError = smi_core.ChannelRegistration.estimateRegistrationError(...
    RegistrationTransform, Coords1, Coords2);

% Create an interpolant function handle for our registration error.
InterpolantFunction = scatteredInterpolant(...
    Coords1(:, 1), Coords1(:, 2), double(SquaredError), ...
    'natural', 'linear');

% Estimate the squared error on a grid.
[XGrid, YGrid] = meshgrid(...
    (FOV(1)-GridSpacing):GridSpacing:(FOV(3)+GridSpacing), ...
    (FOV(2)-GridSpacing):GridSpacing:(FOV(4)+GridSpacing));
SquaredErrorGrid = InterpolantFunction(XGrid(:), YGrid(:));

% Display the grid of estimated squared errors.
SquaredErrorGrid = reshape(SquaredErrorGrid, size(XGrid));
ImshowCorrection = [-0.5, 0.5];
imshow(SquaredErrorGrid, [], 'Parent', PlotAxes, ...
    'XData', FOV([1, 3])+ImshowCorrection, ...
    'YData', FOV([2, 4])+ImshowCorrection)
colormap(PlotAxes, 'parula')
title(PlotAxes, ...
    'Squared registration error estimated in full field of view')
xlabel(PlotAxes, 'X (pixels)')
ylabel(PlotAxes, 'Y (pixels)')
ColorBar = colorbar(PlotAxes);
ColorBar.Label.String = 'Squared error (Pixels^2)';

% Plot the true points as a reference.
hold(PlotAxes, 'on');
line(PlotAxes, Coords1(:, 1), Coords1(:, 2), ...
    'Color', 'k', 'Marker', '.', 'LineStyle', 'none')


end