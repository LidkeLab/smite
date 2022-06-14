function [PlotAxes] = plotImageStack(PlotAxes, ...
    ImageStack, StackSpacing, ColorMap, ViewLOS, ImageScaling)
%plotImageStack plots a stack of images.
% This method will plot a stack of images contained in a 3D array. The
% intention is that these plots will be used for qualitative purposes 
% (e.g., displaying an example of a dataset).
%
% INPUTS:
%   PlotAxes: Axes in which the image stack will be plotted.
%             (Default = gca())
%   ImageStack: A 3d array containing the images. 
%               (MxNx3xNImages numeric array)
%   StackSpacing: Spacing between the NImages in ImageStack.
%                 (pixels)(Default = 5)
%   ColorMap: Color map used to display the images. Note that this is not
%             defined/used the same way as colormaps are sometimes used in
%             MATLAB, and instead defines a splitting ratio between the RGB
%             channels. The use of a fourth element will allow for setting
%             the 'FaceAlpha' property of the surface() call used below.
%             (1x3(4) or NImagesx3(4) numeric array)(Default = [1, 1, 1])
%   ViewLOS: Line of site for viewing the stack. This will be passed
%            directly to the MATLAB method view(). (Default = [-45, 10])
%   ImageScaling: Scaling mode of the images.
%                 'global': Each image is rescaled with respect to the
%                           entire stack.
%                 'local': Each image is independently rescaled to [0, 1].
%                 (char array)(Default = 'global')
%
% OUTPUTS:
%   PlotAxes: Axes in which the image stack was plotted.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Set defaults if needed.
if (~exist('PlotAxes', 'var') || isempty(PlotAxes) || ~isvalid(PlotAxes))
    PlotAxes = gca();
end
if (~exist('StackSpacing', 'var') || isempty(StackSpacing))
    StackSpacing = 5;
end
if (~exist('ColorMap', 'var') || isempty(ColorMap))
    ColorMap = [1, 1, 1];
end
if (~exist('ViewLOS', 'var') || isempty(ViewLOS))
    ViewLOS = [-45, 10];
end
if (~exist('ImageScaling', 'var') || isempty(ImageScaling))
    ImageScaling = 'global';
end

% Reshape input stack if needed.
StackSize = size(ImageStack, 1:4);
if (StackSize(3) ~= 3)
    if (StackSize(3) ~= 1)
        % Assume images are arranged as YxXxNImages.
        ImageStack = reshape(ImageStack, ...
            [StackSize(1), StackSize(2), 1, StackSize(3)]);
    end
    ImageStack = repmat(ImageStack, [1, 1, 3]);
end
NImages = size(ImageStack, 4);

% Reshape 'ColorMap' if needed.
if (size(ColorMap, 1) ~= NImages)
    ColorMap = repmat(ColorMap(1, 1:size(ColorMap, 2)), NImages, 1);
end
FaceAlpha = ones(NImages, 1);
if (size(ColorMap, 2) > 3)
    FaceAlpha = ColorMap(:, 4);
end

% Prepare the plot axes.
XRange = [1, StackSize(2)];
YRange = [1, StackSize(1)];
ZRange = [1, NImages * StackSpacing];
axis(PlotAxes, [XRange, YRange, ZRange])
view(PlotAxes, ViewLOS)
axis(PlotAxes, 'off')
axis(PlotAxes, 'equal')

% Rescale the image stack.
ImageStack = (ImageStack-min(ImageStack(:))) ...
    ./ max(ImageStack(:)-min(ImageStack(:)));
if strcmpi(ImageScaling, 'local')
    for nn = 1:NImages
        ImageStack(:, :, :, nn) = ...
            smi_vis.contrastStretch(ImageStack(:, :, :, nn));
    end
end

% Loop through the image stack and plot each image.
for nn = 1:NImages
    ZSurfVal = repmat(nn * StackSpacing, [2, 2]);
    surface(PlotAxes, XRange, YRange, ZSurfVal, ...
        ImageStack(:, :, :, nn).*reshape(ColorMap(nn, 1:3), 1, 1, 3), ...
        'facecolor', 'texturemap', 'FaceAlpha', FaceAlpha(nn));
end


end