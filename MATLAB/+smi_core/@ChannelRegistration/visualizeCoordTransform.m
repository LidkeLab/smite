function [PlotFigure] = visualizeCoordTransform(PlotFigure, ...
    RegistrationTransform, FrameSize, GridSpacing)
%visualizeCoordTransform creates visuals for a coordinate transform.
% This method will plot the magnitude and gradient of RegistrationTransform
% at points within a FrameSize(1)xFrameSize(2) grid.
%
% NOTE: Default values can be used (when available) by setting the
%       corresponding input to [].  For example, 
%       visualizeCoordTransform([], TForm, FrameSize) will use a default
%       value for PlotFigure.
%
% INPUTS:
%   PlotFigure: Figure in which we will make the plots. (Default = gcf())
%   RegistrationTransform: RegistrationTransform is a MATLAB tform object
%                          describing the transform to be visualized.
%                          (tform object)
%   FrameSize: Size of the frame (ROI) within which RegistrationTransform 
%              can be applied. (Pixels)(2x1, [YSize, XSize])
%   GridSpacing: Spacing of a grid made across FrameSize.  The grid will
%                consist of points
%                [0:GridSpacing:FrameSize(1), 0:GridSpacing:FrameSize(2)] 
%                at which we will apply RegistrationTransform.
%                (Pixels)(Default = 1)
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
if (~exist('GridSpacing', 'var') || isempty(GridSpacing))
    GridSpacing = 1;
end

% Compute a "heatmap" of the transform magnitude across the ROI, as well as
% the gradient of this heatmap.
[XGrid, YGrid] = meshgrid(0:GridSpacing:FrameSize(2), ...
    0:GridSpacing:FrameSize(1));
TransformedCoords = smi_core.ChannelRegistration.transformCoords(...
    RegistrationTransform, [XGrid(:), YGrid(:)]);
XGridTransformed = reshape(TransformedCoords(:, 1), size(XGrid));
YGridTransformed = reshape(TransformedCoords(:, 2), size(YGrid));
TransformMagnitude = sqrt((XGrid-XGridTransformed).^2 ...
    + (YGrid-YGridTransformed).^2);
[XDerivative, YDerivative] = gradient(TransformMagnitude);
GradientMagnitude = sqrt(XDerivative.^2 + YDerivative.^2);
TransformMean = mean(TransformMagnitude(:));
TransformStD = std(TransformMagnitude(:));
OutlierBool = ((TransformMagnitude-TransformMean) > 2*TransformStD);
TransformMagnitude(OutlierBool) = max(TransformMagnitude(~OutlierBool));

% Plot the heatmap.
PlotAxesTop = subplot(1, 2, 1, 'Parent', PlotFigure);
PlotAxesTop.YDir = 'reverse';
axis(PlotAxesTop, 'tight');
axis(PlotAxesTop, 'equal');
hold(PlotAxesTop, 'on');
surface(PlotAxesTop, XGrid, YGrid, TransformMagnitude, ...
    'EdgeColor', 'none')
title(PlotAxesTop, 'sqrt(XCorrection^2 + YCorrection^2)')
xlabel(PlotAxesTop, 'X (pixels)')
ylabel(PlotAxesTop, 'Y (pixels)')
RegColorBar = colorbar(PlotAxesTop);
RegColorBar.Label.String = 'Magnitude of correction (pixels)';
PlotAxesBottom = subplot(1, 2, 2, 'Parent', PlotFigure);
PlotAxesBottom.YDir = 'reverse';
axis(PlotAxesBottom, 'tight');
axis(PlotAxesBottom, 'equal');
hold(PlotAxesBottom, 'on');
surface(PlotAxesBottom, XGrid, YGrid, GradientMagnitude, ...
    'EdgeColor', 'none')
title(PlotAxesBottom, 'sqrt(XDerivative^2 + YDerivative^2)')
xlabel(PlotAxesBottom, 'X (pixels)')
ylabel(PlotAxesBottom, 'Y (pixels)')
RegColorBar = colorbar(PlotAxesBottom);
RegColorBar.Label.String = 'Magnitude of gradient';


end