function saveResults(obj, SaveParams)
%saveResults saves useful results of diffusion estimation analysis.
% This method can be used to save several pieces of information from the
% diffusion analysis results.
%
% INPUTS:
%   SaveParams: Structure array of various parameters defining what we
%               should save.
%               Fields: MakeFitPlot: Generates and saves a plot of the MSD
%                                    fit. (Default = true)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Define the path to the .mat file in which results will be saved.
if isempty(obj.BaseSaveName)
    obj.BaseSaveName = smi_helpers.genTimeString('_');
end
if ~exist(obj.SaveDir, 'dir')
    mkdir(obj.SaveDir);
end

% Check if the SaveParams structure was provided, initializing it to empty
% if it wasn't.
if (~exist('SaveParams', 'var') || isempty(SaveParams))
    SaveParams = struct();
end

% Set default parameter values where needed.
DefaultSaveParams.MakeFitPlot = true;
DefaultParameterNames = fieldnames(DefaultSaveParams);
InputParameterNames = fieldnames(SaveParams);
for ii = 1:numel(DefaultParameterNames)
    if ~any(ismember(DefaultParameterNames{ii}, InputParameterNames))
        % The field DefaultParameterNames{ii} is not present in the
        % SaveParams structure and so the default must be added.
        SaveParams.(DefaultParameterNames{ii}) = ...
            DefaultSaveParams.(DefaultParameterNames{ii});
    end
end

% Save some class properties to a .mat file.
FilePath = fullfile(obj.SaveDir, ...
    ['DiffusionResults_', obj.BaseSaveName, '.mat']);
if (obj.Verbose > 1)
    fprintf('saveResults(): saving diffusion estimates to \n\t%s...\n', ...
        FilePath)
elseif (obj.Verbose > 0)
    fprintf('saveResults(): saving diffusion estimation results...\n')
end
DiffusionStruct = obj.DiffusionStruct;
FitMethod = obj.FitMethod;
MSDEnsemble = obj.MSDEnsemble;
MSDSingleTraj = obj.MSDSingleTraj;
MaxFrameLag = obj.MaxFrameLag;
save(FilePath, 'DiffusionStruct', 'FitMethod', 'MaxFrameLag', ...
    'MSDEnsemble', 'MSDSingleTraj', '-v7.3');

% Generate an MSD fit plot.
if SaveParams.MakeFitPlot
    FilePath = fullfile(obj.SaveDir, ...
        ['MSDEnsembleFit_', obj.BaseSaveName, '.mat']);
    if (obj.Verbose > 1)
        fprintf('saveResults(): saving MSD fit plot to \n\t%s...\n', ...
            FilePath)
    end
    VisiblePlot = smi_helpers.stringMUX({'off', 'on'}, obj.Verbose);
    PlotFigure = figure('Visible', VisiblePlot);
    PlotAxes = axes(PlotFigure);
    obj.plotEnsembleMSD(PlotAxes, ...
        obj.MSDEnsemble, obj.DiffusionStruct, obj.UnitFlag);
    saveas(PlotFigure, FilePath, 'png')
    if (obj.Verbose < 2)
        % For the higher verbosity levels, we'll keep the plot open for the
        % user to view.
        close(PlotFigure)
    end
end


end