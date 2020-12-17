function [PlotAxes] = visualizeImageTransform(PlotAxes, ...
    RegistrationTransform, FrameSize, GridSpacing)
%visualizeImageTransform creates visuals for an image transform.
% This method will create a grid image and apply RegistrationTransform to
% it so as to help us visualize what the transform is doing.
%
% NOTE: Default values can be used (when available) by setting the
%       corresponding input to [].  For example, 
%       visualizeImageTransform([], TForm, FrameSize) will use a default
%       value for PlotAxes.
%
% INPUTS:
%   PlotAxes: Axes in which we will make the plots. (Default = gca())
%   RegistrationTransform: RegistrationTransform is a MATLAB tform object
%                          describing the transform to be visualized.
%                          (tform object)
%   FrameSize: Size of the frame (ROI) within which RegistrationTransform 
%              can be applied. (Pixels)(1x2, [YSize, XSize])
%   GridSpacing: Spacing between grid lines in the grid image that we will
%                transform.  For example, in 1D, GridSpacing = 2 would give
%                us a grid image that looks like [1, 0, 1, 0, ...]
%                (Pixels)(Default = 10)
%
% OUTPUTS:
%   PlotAxes: A MATLAB axes handle for the axes in which the 
%             visualizations are contained. 

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Set default parameters and modify inputs if needed.
if (~exist('PlotAxes', 'var') || isempty(PlotAxes) ...
        || ~isvalid(PlotAxes))
    PlotAxes = gca();
end
if (~exist('GridSpacing', 'var') || isempty(GridSpacing))
    GridSpacing = 10;
end
if (size(FrameSize,  1) > size(FrameSize, 2))
    FrameSize = FrameSize.';
end

% Create a grid image that we can transform.
GridImage = zeros(FrameSize);
GridImage(1:GridSpacing:FrameSize(1), :) = 1;
GridImage(:, 1:GridSpacing:FrameSize(2)) = 1;

% Transform the grid image.
TransformedGrid = smi_core.ChannelRegistration.transformImages(...
    RegistrationTransform, GridImage);

% Create an RGB image containing the original and transformed grid (to help
% show the changes).
ComparisonImage  = zeros([FrameSize, 3]);
ComparisonImage(:, :, 2) = GridImage;
ComparisonImage(:, :, 1) = TransformedGrid;
ComparisonImage(:, :, 3) = TransformedGrid;

% Display the grid image in the desired PlotAxes.
imshow(ComparisonImage, 'Parent', PlotAxes)


end