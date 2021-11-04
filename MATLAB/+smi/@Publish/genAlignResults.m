function [AlignResultsStruct] = genAlignResults(obj, FilePath, SaveDir)
%genAlignResults produces figures/info related to brightfield registration.
% Create figures associated with the initial brightfield registration of a
% cell process on the sequential microscope.
%
% INPUTS:
%   FilePath: Full file path to the location of the .h5 file containing the
%             raw data.
%   SaveDir: Directory associated with Dataset in which the analysis
%            results will be saved.
%
% OUTPUTS:
%   AlignResultsStruct: A structured array containing the alignment results
%                       computed in this method.
%
% REQUIRES:
%   MATLAB Image Processing Toolbox 2014a or later.

% Created by:
%   David J. Schodt (Lidke Lab, 2018)


% Pre-allocate the AlignResultsStruct structured array.
AlignResultsStruct = struct();

% Grab the alignment registration structure from the dataset.
AlignReg = smi_core.LoadData.readH5File(FilePath, 'AlignReg');

% Load the SMR structure for this dataset.
FileNameStruct = dir(fullfile(SaveDir, '*Results.mat'));
if (~isempty(FileNameStruct) ...
        && isfile(fullfile(FileNameStruct.folder, FileNameStruct.name)))
    load(fullfile(SaveDir, FileNameStruct.name), 'SMD');
else
    SMD = struct([]);
end

% If the AlignReg structure is empty, do not proceed.
if isempty(AlignReg)
    return
end

% If needed, make the requested directory.
if ~isfolder(SaveDir)
    mkdir(SaveDir)
end

% Grab the attributes and data from the AlignReg structure, assuming the
% attributes of the first dataset the same for all datasets.
AlignRegData = {AlignReg.Data};

% Grab the error signal history from the last dataset.
ErrorSignalHistory = AlignRegData{end}.ErrorSignalHistory;


% Plot the cumulative error signal during the acquisition, skipping the
% first correction since it is expected to be quite large.
FigureHandle = figure();
PlotAxes = axes(FigureHandle);
NCorrections = size(ErrorSignalHistory, 1);
XPlotArray = (2:NCorrections).';
ErrorSignalSum = cumsum(ErrorSignalHistory(2:end, :), 1);
stairs(PlotAxes, XPlotArray, ErrorSignalSum(:, 1)*1e3)
hold(PlotAxes, 'on');
axis(PlotAxes, 'tight');
stairs(PlotAxes, XPlotArray, ErrorSignalSum(:, 2)*1e3)
stairs(PlotAxes, XPlotArray, ErrorSignalSum(:, 3)*1e3)
plot(PlotAxes, PlotAxes.XLim, [0, 0], 'k:')
title(PlotAxes, 'Cumulative correction')
xlabel(PlotAxes, 'Correction Number')
ylabel(PlotAxes, 'Cumulative error Signal (nm)')
legend(PlotAxes, {'X', 'Y', 'Z', '0 nm reference'}, 'Location', 'best')
saveas(FigureHandle, fullfile(SaveDir, 'AlignRegHistorySum.png'), 'png');
close(FigureHandle);

% Plot the cumulative correction for each dataset.
FigureHandle = figure();
PlotAxes = axes(FigureHandle);
NDatasets = numel(AlignRegData);
Correction = NaN(NDatasets, 3);
Correction(1, :) = sum(AlignRegData{1}.ErrorSignalHistory, 1);
line(PlotAxes, 1, Correction(1, 1)*1e3, ...
    'Marker', 'x', 'Color', [0 0.4470 0.7410], 'LineWidth', 2)
line(PlotAxes, 1, Correction(1, 2)*1e3, ...
    'Marker', 'o', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2)
line(PlotAxes, 1, Correction(1, 3)*1e3, ...
    'Marker', '*', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2)
for ii = 2:NDatasets
    Correction(ii, :) = sum(AlignRegData{ii}.ErrorSignalHistory, 1) ...
        - sum(Correction(1:(ii-1), :), 1);
    line(PlotAxes, ii, Correction(ii, 1)*1e3, ...
        'Marker', 'x', 'Color', [0 0.4470 0.7410], 'LineWidth', 2)
    line(PlotAxes, ii, Correction(ii, 2)*1e3, ...
        'Marker', '*', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2)
    line(PlotAxes, ii, Correction(ii, 3)*1e3, ...
        'Marker', 'o', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2)
end
DummyLines(1) = line(PlotAxes, NaN, NaN, ...
    'LineStyle', 'none', 'Marker', 'x', 'Color', [0 0.4470 0.7410]);
DummyLines(2) = line(PlotAxes, NaN, NaN, ...
    'LineStyle', 'none', 'Marker', '*', 'Color', [0.8500 0.3250 0.0980]);
DummyLines(3) = line(PlotAxes, NaN, NaN, ...
    'LineStyle', 'none', 'Marker', 'o', 'Color', [0.9290 0.6940 0.1250]);
title(PlotAxes, 'Total correction per dataset')
xlabel(PlotAxes, 'Dataset number')
ylabel(PlotAxes, 'Error Signal (nm)')
legend(PlotAxes, DummyLines, {'X', 'Y', 'Z'}, ...
    'Location', 'best')
saveas(FigureHandle, fullfile(SaveDir, 'AlignRegErrorPerDataset.png'), 'png');
close(FigureHandle);

% Plot the error signal during the acquisition, skipping the first
% correction since it is expected to be quite large.
FigureHandle = figure();
PlotAxes = axes(FigureHandle);
NCorrections = size(ErrorSignalHistory, 1);
XPlotArray = (2:NCorrections).';
ErrorSignalPlotArray = ErrorSignalHistory(2:end, :);
plot(PlotAxes, XPlotArray, ErrorSignalPlotArray(:, 1)*1e3)
hold(PlotAxes, 'on');
axis(PlotAxes, 'tight');
plot(PlotAxes, XPlotArray, ErrorSignalPlotArray(:, 2)*1e3)
plot(PlotAxes, XPlotArray, ErrorSignalPlotArray(:, 3)*1e3)
plot(PlotAxes, PlotAxes.XLim, [0, 0], 'k:')
title(PlotAxes, 'Correction history')
xlabel(PlotAxes, 'Correction Number')
ylabel(PlotAxes, 'Error Signal (nm)')
legend(PlotAxes, {'X', 'Y', 'Z', '0 nm reference'}, 'Location', 'best')
saveas(FigureHandle, fullfile(SaveDir, 'AlignRegErrorSignal.png'), 'png');
close(FigureHandle);

% Create interesting movies from the AlignReg data.
[ImagesStruct] = obj.genAlignMovies(AlignRegData, SaveDir);

% Create plots about statistics/information that is derived from the
% alignment result images.
[StatsStruct] = obj.genAlignStats(AlignReg, SMD, SaveDir);

% Create interesting plots from the AlignReg data related to the
% cross-correlation process.
[XCorrStruct] = obj.genAlignXCorr(AlignReg, SaveDir);

% Populate the AlignResultsStruct.
AlignResultsStruct.DiffImages = ImagesStruct.DiffImages;
AlignResultsStruct.OverlayImages = ImagesStruct.OverlayImages;
AlignResultsStruct.SSIM = StatsStruct.SSIM;
AlignResultsStruct.RegError = StatsStruct.RegError;
AlignResultsStruct.MaxCorr = XCorrStruct.MaxCorr;
AlignResultsStruct.MaxCorrFit = XCorrStruct.MaxCorrFit;
AlignResultsStruct.OffsetFitSuccess = XCorrStruct.OffsetFitSuccess;
AlignResultsStruct.MaxIterReached = XCorrStruct.MaxIterReached;


end