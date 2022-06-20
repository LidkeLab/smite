function [ImagesStruct] = genAlignMovies(AlignRegData, SaveDir)
%genAlignMovies generates movies related to brightfield registration.
% This method will create movies which help to visualize the brightfield
% alignment procedure of an acquisition.  Two movies are produced by this
% method.  One is a movie showing overlay images of the in focus
% brightfieldimage reference and the brightfield image at the chosen focus
% after brightfield alignment, and the other contains the scaled difference
% images between the two previously mentioned images.
%
% INPUTS:
%   AlignRegData: Cell array containing the data structures for each
%                 dataset.
%   SaveDir: Directory associated with Dataset in which the analysis
%            results will be saved.
%
% OUTPUTS:
%   ImagesStruct: A structured array containing difference images and
%                 overlay images computed in this method.

% Created by:
%   David J. Schodt (Lidke Lab, 2018)


% Pre-allocate the ImagesStruct;
ImagesStruct = struct();

% Grab the reference image from the last data structure.
ReferenceImage = AlignRegData{end}.Image_Reference;

% Create a stack of the images taken after each registration step has
% completed (i.e. the final brightfield images taken just before proceeding
% with a fluorescence sequence).
NDatasets = numel(AlignRegData);
StackSize = [size(ReferenceImage), NDatasets];
CurrentImages = zeros(StackSize);
for ii = 1:NDatasets
    CurrentImages(:, :, ii) = AlignRegData{ii}.Image_Current;
end

% Prepare video writer objects and set display parameters.
islinux = isunix && ~ismac;
if islinux
   VideoObjectOverlays = VideoWriter(...
       fullfile(SaveDir, 'AlignRegOverlayMovie.avi'), ...
       'Motion JPEG AVI');
   VideoObjectDiffImages = VideoWriter(...
       fullfile(SaveDir, 'AlignRegDiffImageMovie.avi'), ...
       'Motion JPEG AVI');
else
   VideoObjectOverlays = VideoWriter(...
       fullfile(SaveDir, 'AlignRegOverlayMovie.mp4'), ...
       'mpeg-4');
   VideoObjectDiffImages = VideoWriter(...
       fullfile(SaveDir, 'AlignRegDiffImageMovie.mp4'), ...
       'mpeg-4');
end
VideoObjectOverlays.FrameRate = 1; % fps
VideoObjectDiffImages.FrameRate = 1;
VideoObjectOverlays.Quality = 100;
VideoObjectDiffImages.Quality = 100;
VideoObjectOverlays.open();
VideoObjectDiffImages.open();

% Create the movie of the brightfield overlay images after the
% registration corrections.
OverlayImages = zeros([StackSize(1:2), 3, StackSize(3)], 'uint8');
FigureHandle = figure();
PlotAxes = axes(FigureHandle);
for ii = 1:NDatasets
    ImageCurrentFrame = CurrentImages(:, :, ii);
    OverlayImages(:, :, :, ii) = imfuse(...
        ReferenceImage, ImageCurrentFrame, ...
        'falsecolor', 'Scaling', 'independent', ...
        'ColorChannels', [1, 2, 0]);
    imshow(OverlayImages(:, :, :, ii), [], 'Parent', PlotAxes);
    FrameNumDisplay = text(PlotAxes, 10, 10, num2str(ii));
    FrameNumDisplay.Color = [1, 1, 1];
    Frame = getframe(PlotAxes);
    VideoObjectOverlays.writeVideo(Frame);
end
VideoObjectOverlays.close()
cla(PlotAxes)

% Create the movie of the difference images (the difference between the
% scaled reference image and the scaled post correction image).
DiffImages = zeros(StackSize);
for ii = 1:NDatasets
    ImageCurrentFrame = CurrentImages(:, :, ii);
    DiffImages(:, :, ii) = imfuse(ReferenceImage, ImageCurrentFrame, ...
        'diff', 'Scaling', 'independent');
    imshow(DiffImages(:, :, ii), [], 'Parent', PlotAxes);
    FrameNumDisplay = text(PlotAxes, 10, 10, num2str(ii));
    FrameNumDisplay.Color = [1, 1, 1];
    Frame = getframe(PlotAxes);
    VideoObjectDiffImages.writeVideo(Frame);
end
VideoObjectDiffImages.close()
close(FigureHandle)

% Save the overlay and difference image arrays as fields in the
% ImagesStruct.
ImagesStruct.DiffImages = DiffImages;
ImagesStruct.OverlayImages = OverlayImages;


end
