function [] = makeOverlayPlots(ImageShift, RegError, MaxCorr, ...
    SRPixelSize, BPPixelSize, SaveDir)
%makeOverlayPlots makes interesting plots from two color overlays
% This method will estimate shifts and registration errors from the data
% used to generate the two images in ImageStack.
%
% INPUTS:
%   ImageShift: A two column array with each row containing the [row, col]
%               offsets between two images as determined from findshift().
%               (SR pixels)
%   RegError: The estimated brightfield registration errors for datasets
%             used to generate the images corresponding to the shift
%             contained in ImageShift.  RegError is a cell array, with the
%             row index of the cell array corresponding to a row in
%             ImageShift and the column index corresponding to the label
%             number. (raw data pixels)
%   MaxCorr: The maximum correlation coefficients from brightfield
%            registration for the datasets used to generate the second
%            label in the overlay images in consideration.  MaxCorr is a
%            cell array, with the row index of the cell array corresponding
%            to a row in ImageShift and the column index corresponding to
%            the label number.
%   SRPixelSize: Size of a pixel in the SR images as projected onto the
%                sample plane. (microns)(default = 0.005)
%   BPPixelSize: Size of a pixel in the raw data as projected onto the
%                sample plane. (microns)(default = 0.1)
%   SaveDir: Directory in which output plots will be saved.

% Created by:
%   David J. Schodt (Lidke Lab, 2018)


% Set defaults if needed.
if (~exist('SRPixelSize', 'var') || isempty(SRPixelSize))
    SRPixelSize = 0.005; % microns
end
if (~exist('BPPixelSize', 'var') || isempty(BPPixelSize))
    BPPixelSize = 0.1; % microns
end

% Convert ImageShift and RegError to units of nanometers.
ImageShift = ImageShift * SRPixelSize * 1e3; % pixel -> micron -> nanometer
RegError = cellfun(@(x) x * BPPixelSize * 1e3, RegError, ...
    'UniformOutput', false); % pixel -> micron -> nanometer

% Make a scatterplot of the image shifts.
FigureHandle = figure();
PlotAxes = axes(FigureHandle);
hold(PlotAxes, 'on');
NOverlays = size(ImageShift, 1); % number of overlay images considered
ColorMap = colormap(PlotAxes, parula(NOverlays));
for ii = 1:NOverlays
    % Plot the shift between the two images.
    plot(PlotAxes, ImageShift(ii, 2), ImageShift(ii, 1), ...
        'Marker', 'x', 'MarkerEdgeColor', ColorMap(ii, :), ...
        'MarkerSize', 10, 'LineWidth', 1, 'LineStyle', 'none')
end
ColorBar = colorbar(PlotAxes); % add a color bar to show image number
ColorBar.Ticks = linspace(min(ColorBar.Limits), max(ColorBar.Limits), ...
    NOverlays);
ColorBar.TickLabels = 1:NOverlays;
ColorBar.Label.String = 'Image Number';
plot(PlotAxes, [0, 0], PlotAxes.YLim, 'k--')
plot(PlotAxes, PlotAxes.XLim, [0, 0], 'k--')
xlabel(PlotAxes, 'X (nm)')
ylabel(PlotAxes, 'Y (nm)')
title(PlotAxes, 'Image Shifts')

% Save the scatterplot of ImageShift and RegError.
saveas(FigureHandle, ...
    fullfile(SaveDir, 'ImageShiftScatterPlot.png'), 'png');
saveas(FigureHandle, ...
    fullfile(SaveDir, 'ImageShiftScatterPlot.fig'));
close(FigureHandle);

% Plot the magnitude of ImageShift vs. MaxCorr as well as the magnitude of
% ImageShift vs. the magnitude of RegError.
if ~(isempty(MaxCorr) || isempty(RegError))
    FigureHandle = figure();
    SubplotAxesTop = subplot(2, 1, 1, 'Parent', FigureHandle);
    hold(SubplotAxesTop, 'on');
    SubplotAxesBottom = subplot(2, 1, 2, 'Parent', FigureHandle);
    hold(SubplotAxesBottom, 'on');
    for ii = 1:NOverlays
        % Compute the magnitude of the current ImageShift and make a repeated
        % array of this value (to plot vs. the longer arrays MaxCorr and
        % RegError).
        CurrentImageShiftMag = sqrt((ImageShift(ii, 1))^2 ...
            + (ImageShift(ii, 2))^2);
        CurrentImageShiftMag = repelem(CurrentImageShiftMag, ...
            numel(MaxCorr{ii, 1})).';
        
        % Plot the magnitude of ImageShift vs. the MaxCorr.
        plot(SubplotAxesTop, CurrentImageShiftMag, MaxCorr{ii, 1}, ...
            'Marker', 'x', 'MarkerEdgeColor', ColorMap(ii, :), ...
            'LineStyle', 'none')
        plot(SubplotAxesTop, CurrentImageShiftMag, MaxCorr{ii, 2}, ...
            'Marker', 'o', 'MarkerEdgeColor', ColorMap(ii, :), ...
            'LineStyle', 'none')
        
        % Compute the magnitude of the current RegError.
        CurrentRegErrorLabel1 = RegError{ii, 1};
        CurrentRegErrorMagLabel1 = sqrt((CurrentRegErrorLabel1(:, 1)).^2 ...
            + (CurrentRegErrorLabel1(:, 2)).^2);
        CurrentRegErrorLabel2 = RegError{ii, 2};
        CurrentRegErrorMagLabel2 = sqrt((CurrentRegErrorLabel2(:, 1)).^2 ...
            + (CurrentRegErrorLabel2(:, 2)).^2);
        
        % Plot the magnitude of RegError vs. the magnitude of ImageShift.
        plot(SubplotAxesBottom, ...
            CurrentImageShiftMag, CurrentRegErrorMagLabel1, ...
            'Marker', 'x', 'MarkerEdgeColor', ColorMap(ii, :), ...
            'LineStyle', 'none')
        plot(SubplotAxesBottom, ...
            CurrentImageShiftMag, CurrentRegErrorMagLabel2, ...
            'Marker', 'o', 'MarkerEdgeColor', ColorMap(ii, :), ...
            'LineStyle', 'none')
    end
    ColorBar = colorbar(SubplotAxesTop); % add a color bar to show image number
    ColorBar.Ticks = linspace(min(ColorBar.Limits), max(ColorBar.Limits), ...
        NOverlays);
    ColorBar.TickLabels = 1:NOverlays;
    ColorBar.Label.String = 'Image Number';
    ColorBar = colorbar(SubplotAxesBottom); % add a color bar
    ColorBar.Ticks = linspace(min(ColorBar.Limits), max(ColorBar.Limits), ...
        NOverlays);
    ColorBar.TickLabels = 1:NOverlays;
    ColorBar.Label.String = 'Image Number';
    CrossPlot = plot(SubplotAxesTop, NaN, NaN, 'kx'); % 'dummy' points
    CirclePlot = plot(SubplotAxesTop, NaN, NaN, 'ko');
    legend(SubplotAxesTop, [CrossPlot, CirclePlot], {'Label 1', 'Label 2'}, ...
        'Location', 'best');
    CrossPlot = plot(SubplotAxesBottom, NaN, NaN, 'kx'); % 'dummy' points
    CirclePlot = plot(SubplotAxesBottom, NaN, NaN, 'ko');
    legend(SubplotAxesBottom, [CrossPlot, CirclePlot], {'Label 1', 'Label 2'}, ...
        'Location', 'best');
    xlabel(SubplotAxesTop, 'Image Shift Magnitude (nm)')
    xlabel(SubplotAxesBottom, 'Image Shift Magnitude (nm)')
    ylabel(SubplotAxesTop, 'Max. Corr. Coeff.')
    ylabel(SubplotAxesBottom, 'Reg. Error Magnitude (nm)')
    
    % Save the plots of ImageShift vs. MaxCorr/RegError.
    saveas(FigureHandle, ...
        fullfile(SaveDir, 'ImageShiftRegErrorMaxCorr.png'), 'png');
    saveas(FigureHandle, ...
        fullfile(SaveDir, 'ImageShiftRegErrorMaxCorr.fig'));
    close(FigureHandle);
    
    % Plot the magnitude of RegError vs. the cross-correlation maxima.
    FigureHandle = figure();
    PlotAxes = axes(FigureHandle);
    hold(PlotAxes, 'on');
    for ii = 1:NOverlays
        CurrentRegErrorLabel1 = RegError{ii, 1};
        CurrentRegErrorMagLabel1 = sqrt((CurrentRegErrorLabel1(:, 1)).^2 ...
            + (CurrentRegErrorLabel1(:, 2)).^2);
        CurrentRegErrorLabel2 = RegError{ii, 2};
        CurrentRegErrorMagLabel2 = sqrt((CurrentRegErrorLabel2(:, 1)).^2 ...
            + (CurrentRegErrorLabel2(:, 2)).^2);
        plot(PlotAxes, MaxCorr{ii, 1}, CurrentRegErrorMagLabel1, 'bx')
        plot(PlotAxes, MaxCorr{ii, 2}, CurrentRegErrorMagLabel2, 'bo')
    end
    CrossPlot = plot(PlotAxes, NaN, NaN, 'bx'); % 'dummy' points for legend
    CirclePlot = plot(PlotAxes, NaN, NaN, 'bo');
    legend(PlotAxes, [CrossPlot, CirclePlot], {'Label 1', 'Label 2'}, ...
        'Location', 'best');
    xlabel(PlotAxes, 'Cross-Correlation Maximum')
    ylabel(PlotAxes, 'Registration Error Magnitude (nm)')
    title(PlotAxes, 'Registration Errors vs. Cross-Correlation Maxima')
    
    % Save the registration error vs. cross-corr. max.
    saveas(FigureHandle, ...
        fullfile(SaveDir, 'RegErrorVsCrossCorr.png'), 'png');
    saveas(FigureHandle, ...
        fullfile(SaveDir, 'RegErrorVsCrossCorr.fig'));
    close(FigureHandle);
end


end