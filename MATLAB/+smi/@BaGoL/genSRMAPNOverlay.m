function OverlayImageCircle ...
  = genSRMAPNOverlay(SMD, MAPN, XSize, YSize, PixelSize, ...
                      SaveDir, Xstart, Ystart, RadiusScale, ScaleBarLength)
%genSRMAPNOverlay generates a multicolor overlay containing circles with radii
% proportional to the average localization precision (the generalized mean
% between the X and Y precisions) for each localization from an SMD and a BaGoL
% MAPN structure.  The overlay is saved as a png-file with the identifier
% Overlay_SR_Map_circle.
%
% INPUTS:
%    SMD            SMD structure with the following fields (nm):
%                      X, Y, X_SE, Y_SE
%    MAPN           MAPN structure with the following fields (nm):
%                      X, Y, X_SE, Y_SE
%    XSize          x-size of image (nm)
%    YSize          y-size of image (nm)
%    PixelSize      either pixel linear dimension (nm), or string 'rescale'
%                   where an appropriate pixel size will be calculated.
%    SaveDir        directory in which to save the overlay image file
%    XStart:        X starting coordinate of the output image
%    YStart:        Y starting coordinate of the output image
%    RadiusScale    [OPTIONAL, Default = 1] scalar to multiply the localization
%                   precision by
%    ScaleBarLength [OPTIONAL, Default = 1000] length of scale bar on generated
%                   image (nm)
%
% OUTPUT:
%    OverlayCircleImage   overlay image with scalebar
%                   If genSRMAPNOverlay is called with no output argument,then
%                   this image will be saved in SaveDir as
%                   Overlay_SR_Map_circle.png

% Created by:
%    Michael J. Wester, Mohamadreza Fazel (2021) and David J. Schodt (Lidke Lab, 2018)

if ~exist('RadiusScale', 'var')
   RadiusScale = 1;
end
if ~exist('ScaleBarLength', 'var')
   ScaleBarLength = 1000;
end
if ~exist('Xstart', 'var')
   Xstart = 0;
end
if ~exist('Ystart', 'var')
   Ystart = 0;
end
if ~exist('PixelSize','var') || strcmp(PixelSize,'rescale')
    %If not given, the pixel size will be calculated here based on the
    %smallest precision.
    MinPixelsPerCircle = 16;
    CircleRadius = sqrt((MAPN.X_SE.^2+MAPN.Y_SE.^2) / 2) .* RadiusScale;
    SmallestCircumference= 2 * pi * min(CircleRadius);
    if SmallestCircumference == 0
       error('Precision cannot be zero'); 
    end
    PixelSize = SmallestCircumference/MinPixelsPerCircle;
end

% Compute parameters needed for the circle images (these are computed for
% the MAPN dataset and then left the same for the SMD).
BitDepth = 8; % specific to the save type: this is for a 8 bit png
ImageSize = ceil([XSize, YSize]/PixelSize);
if ((ImageSize(1)*ImageSize(2)*3) > (2^32 - 1))
    % imwrite() won't work for an image this large so we'll
    % have to settle with something smaller.
    % NOTE: The 3 corresponds to the 3 color channels that will
    %       be in the final overlay image.
    ImageSize = floor(2^16 * [1, 1] / sqrt(3));
    PixelSize = max([XSize,YSize])/size(ImageSize,1);
end
fprintf('ImageSize = %d x %d\n', ImageSize);

for ii = 1 : 2
   CircleImage = zeros(ImageSize, 'uint8'); % create the 0 image
   if ii == 1
      X = ((SMD.X-Xstart) / PixelSize);
      Y = ((SMD.Y-Ystart) / PixelSize);
      X_SE = (SMD.X_SE / PixelSize);
      Y_SE = (SMD.Y_SE / PixelSize);
   else
      X = ((MAPN.X-Xstart) / PixelSize);
      Y = ((MAPN.Y-Ystart) / PixelSize);
      X_SE = (MAPN.X_SE / PixelSize);
      Y_SE = (MAPN.Y_SE / PixelSize);
   end

   % Loop through each localization and construct the circle image.
   for mm = 1:numel(X)
       % Create the binary circle image by creating ~2*circumference
       % points for each circle (the extra points help make the circle
       % appear smooth).
       CircleRadius = sqrt((X_SE(mm)^2 + Y_SE(mm)^2) / 2) * RadiusScale;
       Theta = linspace(0, 2*pi, ceil(4 * pi * CircleRadius));
       CircleX = CircleRadius*cos(Theta) + X(mm);
       CircleY = CircleRadius*sin(Theta) + Y(mm);
       for nn = 1:numel(Theta)
           % Place the computed points for each circle in the binary image.
           ImageRow = round(CircleY(nn)); % row in the image
           ImageCol = round(CircleX(nn)); % column in the image
           if ((ImageRow<ImageSize(1)) && (ImageRow>1) ...
                   && (ImageCol<ImageSize(2)) && (ImageCol>1))
               % This point lies within the bounds of the image.
               CircleImage(ImageRow, ImageCol) = 2^BitDepth - 1;
               if ii == 2
                   CircleImages(ImageRow, ImageCol,1) = 0;
               end
           end
       end
   end
   CircleImages(1:ImageSize(1), 1:ImageSize(2), ii) = CircleImage;
end

% Produce a green (G) and magenta (R+B) color overlay for the two
% images, with the first image being green.
% If the input stack wasn't a float, make the output OverlayImage a uint8.
OverlayImageCircle = ...
    zeros([size(CircleImages, 1), size(CircleImages, 2), 3], 'uint8');
% Scale the color channels individually to [0, 1].
OverlayImageCircle(:, :, 1) = CircleImages(:, :, 2) ...
    ./ max(max(CircleImages(:, :, 2))); % red
OverlayImageCircle(:, :, 2) = CircleImages(:, :, 1) ...
    ./ max(max(CircleImages(:, :, 1))); % green
OverlayImageCircle(:, :, 3) = CircleImages(:, :, 2) ...
    ./ max(max(CircleImages(:, :, 2))); % blue
% For non-float images, we need to rescale the pixel values.
if ~isfloat(OverlayImageCircle)
    OverlayImageCircle = OverlayImageCircle * 255;
end

OverlayImageCircleName = sprintf('Overlay_SR_Map_circle.png');
% Add in a scalebar accounting for the size
% AdjustedScaleBarLength = ScaleBarLength * ImageSize(1) / XSize;
OverlayImageCircle = ...
   BaGoL.scalebar(OverlayImageCircle, PixelSize, ScaleBarLength);
if nargout == 0
   imwrite(OverlayImageCircle, jet, fullfile(SaveDir, OverlayImageCircleName));
end

end
