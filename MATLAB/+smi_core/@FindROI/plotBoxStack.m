function Data = plotBoxStack(SMD, Data, BoxSize, Params)
%plotBoxStack plots the found boxes in Data/SMD structure.
% This method extends plotBox() to allow slices through multiple frames by
% interactive use of a slide bar.
%
% INPUTS
%    SMD: Single Molecule Data structure with fields
%         FrameNum: Frame numbers of localizations.
%         XBoxCorner: X coordinates of top right box corners. (pixels)
%         YBoxCorner: Y coordinates of top right box corners. (pixels)
%    Data: Raw data images. (XSize x YSize x NFrames)
%    BoxSize: Linear box size used for fitting. (pixels)
%    Params: Set of parameters that can be used to scale the data. 
%            (see image related parameters in
%            smi_vis.GenerateMovies.prepDefaults() as well as usage below). 
%            (Default = smi_vis.GenerateMovies.prepDefaults())
% OUTPUT
%    Data: Scaled data images

% Created by:
%    David J. Schodt (Lidke Lab, 2022)


% Set defaults.
if (~exist('Params', 'var') || isempty(Params))
    Params = smi_vis.GenerateMovies.prepDefaults();
else
    Params = smi_helpers.padStruct(Params, ...
        smi_vis.GenerateMovies.prepDefaults());
end

% Prepare the data.
Data = smi_vis.contrastStretch(single(Data), [0; 1], ...
    Params.PercentileCeiling, Params.PercentileFloor, ...
    Params.MinScaleIntensity);

% Prepare the slide bar.
PlotFigure = figure();
PlotFigure.Name = 'Raw data ROIs from SMD';
SlidePanel = uipanel(PlotFigure, ...
    'Units', 'normalized', 'Position', [0, 0, 1, 0.1]);
NFrames = size(Data, 3);
uicontrol('Parent', SlidePanel, ...
    'Style', 'slider', ...
    'Units', 'normalized', 'Position', [0, 0.5, 1, 0.5], ...
    'HorizontalAlignment', 'left', ...
    'Min', 1, 'Max', NFrames, ...
    'SliderStep', [1, 10] / NFrames, ...
    'Value', 1, ...
    'Callback', @frameSlider);
FrameNumDisplay = uicontrol('Parent', SlidePanel, ...
    'Style', 'text', 'String', sprintf('Frame %i of %i', 1, NFrames),...
    'Units', 'normalized', 'Position', [0, 0, 1, 0.5], ...
    'HorizontalAlignment', 'left');

% Prepare the data panel.
DataPanel = uipanel(PlotFigure, ...
    'Units', 'normalized', 'Position', [0, 0.1, 1, 0.9]);
PlotAxes = axes(DataPanel);
imshow(Data(:, :, 1), [0, 1], 'Parent', PlotAxes)
hold(PlotAxes, 'on')
axis(PlotAxes, 'tight')

% Isolate some bits of data for speed.
XCorners = double(SMD.XBoxCorner);
YCorners = double(SMD.YBoxCorner);
FrameNum = double(SMD.FrameNum);

% Plot the initial boxes.
OnesArray = ones(BoxSize, 1, 'double');
IndArray = (0:(BoxSize-1)).';
plotBoxes(XCorners, YCorners, BoxSize, FrameNum, 1)

% Define the action of the slide bar.
    function frameSlider(Source, ~)
        % Plot the raw data.
        cla(PlotAxes)
        FrameNumber = round(Source.Value);
        imshow(Data(:, :, FrameNumber), [0, 1], 'Parent', PlotAxes)
        FrameNumDisplay.String = sprintf('Frame %i of %i', FrameNumber, NFrames);

        % Plot the boxes.
        plotBoxes(XCorners, YCorners, BoxSize, FrameNum, FrameNumber)
    end

    function plotBoxes(XCorners, YCorners, BoxSize, FrameNum, FrameNumber)
        % Plot the boxes for the current frame.
        CurrentFrame = find(FrameNum == FrameNumber);
        for nn = CurrentFrame.'
            plot(PlotAxes, XCorners(nn)*OnesArray, ...
                YCorners(nn) + IndArray, 'y')
            plot(PlotAxes, (XCorners(nn)+BoxSize-1)*OnesArray, ...
                YCorners(nn) + IndArray, 'y')
            plot(PlotAxes, XCorners(nn) + IndArray, ...
                YCorners(nn)*OnesArray, 'y')
            plot(PlotAxes, XCorners(nn) + IndArray, ...
                (YCorners(nn)+BoxSize-1)*OnesArray, 'y')
        end
    end


end
