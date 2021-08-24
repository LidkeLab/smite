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
title(PlotAxes, 'History of the error signal')
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