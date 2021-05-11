function [OverlayImage, ColorOrderTag] = overlayNImages(ImageStack)
%overlayNImages will fuse N images together into a single overlay image.
% This method will take a stack of N images (an mxnxN array) and overlay
% these images (each being a different color channel) into an output mxn
% image OverlayImage.  N must be between 2 and 4.
%
% NOTE from DJS 9/1/2020: I seem to remember there being a bug with the
%                         color combination math that just happens to not
%                         show up with the colors chosen below, but I can't
%                         remember how it came up/how to fix it at this
%                         time.
%
% INPUTS:
%   ImageStack: (mxnxN) Array of N mxn images which are to be overlain.
%
% OUTPUTS:
%   OverlayImage: (mxnx3) Three color RGB array of the N images from
%                 ImageStack overlain in different color channels.
%   ColorOrderTag: A character array to specify the ordering of the colors.
%                  For example, if ImageStack is mxnx2, OverlayImage will
%                  be a green and magenta overlay of these two images.  The
%                  image ImageStack(:, :, 1) will correspond to green and
%                  ImageStack(:, :, 2) will correspond to magenta, and thus
%                  ColorOrderTag = 'GM'.
%
% Created by:
%   David J. Schodt (Lidke Lab, 2018)


% Ensure that more than one image exists in the stack.  If not, set
% OverlayImage to the input and do not proceed.
if (size(ImageStack, 3) < 2)
    OverlayImage = ImageStack;
    return
end

% Ensure that there aren't more than 4 images in the stack (this code isn't
% written to handle that many).
if (size(ImageStack, 3) > 4)
    error('Input ImageStack cannot contain more than 4 images.');
end

% Generate our color overlay images (3 channel images).
if isfloat(ImageStack)
    % If the input stack was a float, make the output OverlayImage a
    % single.
    OverlayImage = zeros([size(ImageStack, 1), size(ImageStack, 2), 3], ...
        'single');
else
    % If the input stack wasn't a float, make the output OverlayImage a
    % uint8.
    OverlayImage = zeros([size(ImageStack, 1), size(ImageStack, 2), 3], ...
        'uint8');
end
NImages = size(ImageStack, 3); % number of images to overlay
if (NImages == 2)
    % Produce a green (G) and magenta (R+B) color overlay for the two
    % image, with the first image being green.
    
    % Seperate our images to improve code readability.
    Image1 = ImageStack(:, :, 1);
    Image2 = ImageStack(:, :, 2);
    
    % Define the red channel of the overlay image.
    OverlayImage(:, :, 1) = Image2;
    
    % Define the green channel of the overlay image.
    OverlayImage(:, :, 2) = Image1;
    
    % Define the blue channel of the overlay image.
    OverlayImage(:, :, 3) = Image2;
    
    % Scale the color channels individually to [0, 1].
    OverlayImage(:, :, 1) = OverlayImage(:, :, 1) ...
        ./ max(max(OverlayImage(:, :, 1)));
    OverlayImage(:, :, 2) = OverlayImage(:, :, 2) ...
        ./ max(max(OverlayImage(:, :, 2)));
    OverlayImage(:, :, 3) = OverlayImage(:, :, 3) ...
        ./ max(max(OverlayImage(:, :, 3)));
    
    % Save a filename tag to specify the color order for the images
    % (first character -> color for image 1,
    % second character -> color for image 2, ...).
    ColorOrderTag = 'GM';
end

if (NImages == 3)
    % Produce a green (G), magenta (R+B), and cyan (G+B) color overlay for
    % the three images, with the first images being green, second image
    % magenta, and third image cyan.
    
    % Seperate our images to improve code readability.
    Image1 = ImageStack(:, :, 1);
    Image2 = ImageStack(:, :, 2);
    Image3 = ImageStack(:, :, 3);
    
    % Define the red channel of the overlay image.
    OverlayImage(:, :, 1) = Image2;
    
    % Define the green channel of the overlay image.
    OverlayImage(:, :, 2) = Image1 + Image3;
    
    % Define the blue channel of the overlay image.
    OverlayImage(:, :, 3) = Image2 + Image3;
    
    % Scale the color channels individually to [0, 1].
    OverlayImage(:, :, 1) = OverlayImage(:, :, 1) ...
        ./ max(max(OverlayImage(:, :, 1)));
    OverlayImage(:, :, 2) = OverlayImage(:, :, 2) ...
        ./ max(max(OverlayImage(:, :, 2)));
    OverlayImage(:, :, 3) = OverlayImage(:, :, 3) ...
        ./ max(max(OverlayImage(:, :, 3)));
    
    % Save a filename tag to specify the color order for the
    % images (first character -> color for image 1, second
    % character -> color for image 2, ...).
    ColorOrderTag = 'GMC';
end

if (NImages == 4)
    % Produce a green (G), magenta (R+B), cyan (G+B), and yellow (R+G)
    % color overlay for the three images, with the first image being green,
    % second image magenta, third image cyan, and fourth image yellow.
    
    % Seperate our images to improve code readability.
    Image1 = ImageStack(:, :, 1);
    Image2 = ImageStack(:, :, 2);
    Image3 = ImageStack(:, :, 3);
    Image4 = ImageStack(:, :, 4);
    
    % Define the red channel of the overlay image.
    OverlayImage(:, :, 1) = Image2 + Image4;
    
    % Define the green channel of the overlay image.
    OverlayImage(:, :, 2) = Image1 + Image3 + Image4;
    
    % Define the blue channel of the overlay image.
    OverlayImage(:, :, 3) = Image2 + Image3;
    
    % Scale the color channels individually to [0, 1].
    OverlayImage(:, :, 1) = OverlayImage(:, :, 1) ...
        ./ max(max(OverlayImage(:, :, 1)));
    OverlayImage(:, :, 2) = OverlayImage(:, :, 2) ...
        ./ max(max(OverlayImage(:, :, 2)));
    OverlayImage(:, :, 3) = OverlayImage(:, :, 3) ...
        ./ max(max(OverlayImage(:, :, 3)));
    
    % Save a filename tag to specify the color order for the
    % images (first character -> color for image 1,
    % second character -> color for image 2, ...).
    ColorOrderTag = 'GMCY';
end

% For non-float images, we need to rescale the pixel values.
if ~isfloat(OverlayImage)
    OverlayImage = OverlayImage * 255;
end


end