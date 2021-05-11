function [StatsStruct] = genAlignStats(AlignRegStruct, SMD, SaveDir)
%genAlignStats generates interesting plots from AlignRegStruct.
% This method will create various plot(s) and analysis data related to the
% brightfield registration process of the acquistion.  A plot will be
% created to show information/statistics derived from the difference image
% between the brightfield image reference and the brightfield image taken
% at the chosen focus during the acquisition. This plot will be saved in
% SaveDir.  Another plot will be created to show the structural similarity
% between the full scale histogram stretched reference image and current
% image for each sequence.  This plot will be saved in SaveDir.
%
% INPUTS:
%   AlignRegStruct: Structured array containing the brightfield
%                   registration data structures for each dataset.
%   SMD: Single Molecule Data structure containing data related to the
%        SR analysis.
%   SaveDir: Directory associated with Dataset in which the analysis
%            results will be saved.
%
% OUTPUTS:
%   StatsStruct: A structured array which will contain the computed
%                statistics/interesting metrics for the input
%                AlignRegStruct.
%
% REQUIRES:
%   MATLAB 2018a or later (to use property 'WindowState' of a figure)
%   MATLAB Image Processing Toolbox 2014a or later (to use ssim())

% Created by:
%   David J. Schodt (Lidke Lab, 2018)


% Pre-allocate our StatsStruct structured array.
StatsStruct = struct();

% Determine the number of sequences contained in this dataset.
NSequences = numel(AlignRegStruct);

% Grab the reference image from the data structure (the index shouldn't
% matter, each dataset should have to same reference image saved).
ReferenceImage = AlignRegStruct(1).Data.Image_Reference;

% Perform a full scale histogram stretch on the reference image scale it to
% [0, 1].
ReferenceImage = (ReferenceImage-min(ReferenceImage(:))) ...
    / (max(ReferenceImage(:))-min(ReferenceImage(:)));

% Create the brightfield difference image for each dataset and plot the
% histogram of it's absolute value.
DiffImages = zeros([size(ReferenceImage), NSequences]);
for ii = 1:NSequences
    % Grab the current image for the ii-th sequence.
    CurrentImage = AlignRegStruct(ii).Data.Image_Current;
    
    % Perform a full scale histogram stretch on the current image.
    CurrentImage = (CurrentImage-min(CurrentImage(:))) ...
        / (max(CurrentImage(:))-min(CurrentImage(:)));
    
    % Create the difference images.
    DiffImages(:, :, ii) = (CurrentImage./max(CurrentImage(:))) ...
        - (ReferenceImage./max(ReferenceImage(:)));
end
MinAbsDiff = min(abs(DiffImages(:)));
MaxAbsDiff = max(abs(DiffImages(:)));
PlotFigure = figure();
PlotFigure.WindowState = 'maximized';
for ii = 1:NSequences
    % Create a subplot in the figure for the current histogram.
    PlotAxes = subplot(NSequences, 1, ii, 'Parent', PlotFigure);
    
    % Plot the histogram of pixel values in the difference image.
    DiffImage = DiffImages(:, :, ii);
    StackedDiffImage = DiffImage(:);
    histogram(PlotAxes, abs(StackedDiffImage), 200);
    hold(PlotAxes, 'on');
    
    % Plot a mean line and a standard deviation line of the histogram.
    plot(PlotAxes, [mean(abs(StackedDiffImage)), ...
        mean(abs(StackedDiffImage))], ...
        PlotAxes.YLim, 'r:', 'LineWidth', 2)
    
    % Modify the appearance of the plot.
    PlotAxes.XLim = [MinAbsDiff, MaxAbsDiff];
    legend(PlotAxes, {sprintf('Sequence %i', ii), ...
        sprintf('Mean = %.4f', mean(abs(StackedDiffImage)))}, ...
        'Location', 'best')
    if (ii == 1)
        title(PlotAxes, ...
            ['Histogram of Absolute Pixelwise Difference of ', ...
            'FSHS Images'], 'FontSize', 15)
    end
    if (ii == ceil(NSequences/2))
        ylabel(PlotAxes, 'Counts', 'FontSize', 15)
    end
    if (ii == NSequences)
        xlabel(PlotAxes, 'Pixelwise Difference', 'FontSize', 15)
    end
end

% Save the plot in the SaveDir and then close the figure.
saveas(PlotFigure, fullfile(SaveDir, 'DiffImageHistogram.png'), 'png');
close(PlotFigure);

% Compute the structural similarity between the full scale histogram
% stretched reference and current images.
SSIM = zeros(NSequences, 1);
for ii = 1:NSequences
    % Grab the current image for the ii-th sequence.
    CurrentImage = AlignRegStruct(ii).Data.Image_Current;
    
    % Perform a full scale histogram stretch on the current image.
    CurrentImage = (CurrentImage-min(CurrentImage(:))) ...
        / (max(CurrentImage(:))-min(CurrentImage(:)));
    
    % Compute the SSIM.
    SSIM(ii) = ssim(CurrentImage, ReferenceImage);
end

% Store the SSIM in the StatsStruct.
StatsStruct.SSIM = SSIM;

% Plot the SSIM and save that plot in SaveDir.
PlotFigure = figure();
PlotAxes = axes(PlotFigure);
hold(PlotAxes, 'on')
plot(PlotAxes, 1:NSequences, SSIM, 'x')
plot(PlotAxes, PlotAxes.XLim, ones(2, 1) * mean(SSIM), '--')
legend(PlotAxes, {'SSIM', sprintf('Mean SSIM = %.4f', mean(SSIM))}, ...
    'Location', 'best')
xlabel(PlotAxes, 'Sequence Number')
ylabel(PlotAxes, 'SSIM')
title(PlotAxes, ...
    'Structural Similarity between Reference Image and Current Image')
saveas(PlotFigure, fullfile(SaveDir, 'SSIM.png'), 'png');
close(PlotFigure);

% Compute the 'registration error' as found from the residual drift that is
% corrected during the drift correction process (if the input SMD wasn't
% empty).
if ~isempty(SMD)
    PlotFigure = figure();
    PlotAxes = axes(PlotFigure);
    [~, RegError] = obj.plotXYRegError(PlotAxes, SMD);
    saveas(PlotFigure, fullfile(SaveDir, 'XYRegError.png'), 'png');
    close(PlotFigure)
    
    % Store the registration error array in the StatsStruct.
    StatsStruct.RegError = RegError;
else
    StatsStruct.RegError = [];
end


end