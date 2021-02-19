function saveResults(obj)
%saveResults saves useful results of diffusion estimation analysis.
% This method can be used to save several pieces of information from the
% diffusion analysis results.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)
    

% Define the path to the .mat file in which results will be saved.
if isempty(obj.BaseSaveName)
    obj.BaseSaveName = smi_helpers.genTimeString('_');
end
if ~exist(obj.SaveDir, 'dir')
    mkdir(obj.SaveDir);
end

% Save some class properties to FileName.
DiffusionStruct = obj.DiffusionStruct;
FitMethod = obj.FitMethod;
MSDEnsemble = obj.MSDEnsemble;
MSDSingleTraj = obj.MSDSingleTraj;
MaxFrameLag = obj.MaxFrameLag;
save(fullfile(obj.SaveDir, ['DiffusionResults_', obj.BaseSaveName]), ...
    'DiffusionStruct', 'FitMethod', ...
    'MSDEnsemble', 'MSDSingleTraj', 'MaxFrameLag');

% Generate an MSD fit plot.
VisiblePlot = smi_helpers.stringMUX({'off', 'on'}, obj.Verbose);
PlotFigure = figure('Visible', VisiblePlot);
PlotAxes = axes(PlotFigure);
obj.plotEnsembleMSD(PlotAxes, ...
    obj.MSDEnsemble, obj.DiffusionStruct, obj.UnitFlag);
saveas(PlotFigure, ...
    fullfile(obj.SaveDir, ['MSDEnsembleFit_', obj.BaseSaveName]), 'png')
close(PlotFigure)


end