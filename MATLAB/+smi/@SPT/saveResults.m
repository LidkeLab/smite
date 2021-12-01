function saveResults(obj)
%saveResults saves various tracking results produced by smi.SPT
% This method will save the SMD, TR, and SMF structures in .mat files
% in a Results folder.  This method will also save 2D and 3D plots of the
% tracked trajectories, as well as a 3D movie of the tracked trajectories
% (if requested).

% Created by:
%   David J. Schodt (Lidke Lab, 2018)


% Create the Results directory within the obj.SaveDir (if needed).
if (exist(obj.SMF.Data.ResultsDir, 'dir') ~= 7)
    % exist() will return a 7 if the ResultsDir is an existing directory.
    mkdir(obj.SMF.Data.ResultsDir)
end

% Create the filename for the saved results.
% NOTE: For now, this is only setup for a single file (hence the {1}).
[~, FileName] = fileparts(obj.SMF.Data.FileName{1});
BaseName = [FileName, '_', obj.SMF.Data.AnalysisID];
ResultsFileName = [BaseName, '_Results.mat'];

% Move the data structures of interest into the workspace with appropriate
% names for saving.
SMD = obj.SMD;
SMDPreThresh = obj.SMDPreThresh;
TR = obj.TR;
SMF = obj.SMF;

% Save the data in a .mat file.
save(fullfile(obj.SMF.Data.ResultsDir, ResultsFileName), ...
    'SMD', 'SMDPreThresh', 'TR', 'SMF');

% If pre-channel registration results were stored, we'll save those too.
if ~(isempty(obj.SMDPreCR) || isempty(obj.SMDPreThreshPreCR) ...
        || isempty(obj.TRPreCR))
    ResultsFileNamePreCR = [BaseName, '_Results_PreCR.mat'];
    SMD = obj.SMDPreCR;
    SMDPreThresh = obj.SMDPreThreshPreCR;
    TR = obj.TRPreCR;
    save(fullfile(obj.SMF.Data.ResultsDir, ResultsFileNamePreCR), ...
        'SMD', 'SMDPreThresh', 'TR', 'SMF');
end

% Create a movie of the tracks and save the resulting movie.
if obj.GenerateMovies
    % Load the raw data.
    LD = smi_core.LoadData;
    [~, RawData, obj.SMF] = ...
        LD.loadRawData(obj.SMF, 1, obj.SMF.Data.DataVariable);
    
    % Generate and save the movie.
    if ~(isempty(RawData) || isempty(obj.TR))
        MovieMaker = smi_vis.GenerateMovies(obj.MovieParams);
        MovieMaker.TR = obj.TR;
        MovieMaker.RawData = RawData;
        MovieMaker.SMF = obj.SMF;
        MovieFileName = [BaseName, '_movie.mp4'];
        MovieMaker.saveMovie(fullfile(obj.SMF.Data.ResultsDir, MovieFileName))
    elseif (obj.Verbose > 0)
        warning('smi.SPT.saveResults(): no movie produced, data is empty!')
    end
end

% Create and save 2D and 3D trajectory plots.
if obj.GeneratePlots
    % Make and save the 2D plot.
    if ~isempty(obj.TR)
        PlotFigure = figure();
        PlotAxes = axes(PlotFigure);
        MovieMaker = smi_vis.GenerateMovies(obj.MovieParams);
        MovieMaker.SMF = obj.SMF;
        MovieMaker.TR = obj.TR;
        MovieMaker.setVitalParams()
        MovieMaker.prepAxes(PlotAxes);
        EmptySMD = smi_core.SingleMoleculeData.createSMD();
        MovieMaker.makeFrame(PlotAxes, obj.TR, [], MovieMaker.Params, ...
            obj.SMF, EmptySMD, obj.TR(1).NFrames);
        Traj2DFileName = [BaseName, '_plot2D'];
        saveas(PlotFigure, fullfile(obj.SMF.Data.ResultsDir, Traj2DFileName))
        saveas(PlotFigure, ...
            fullfile(obj.SMF.Data.ResultsDir, [Traj2DFileName, '.png']))
        
        % Make and save the 3D plot.
        MovieMaker.Params.LineOfSite = [-45, 15];
        MovieMaker.prepAxes(PlotAxes);
        EmptySMD = smi_core.SingleMoleculeData.createSMD();
        MovieMaker.makeFrame(PlotAxes, obj.TR, [], MovieMaker.Params, ...
            obj.SMF, EmptySMD, obj.TR(1).NFrames);
        Traj3DFileName = [BaseName, '_plot3D'];
        saveas(PlotFigure, fullfile(obj.SMF.Data.ResultsDir, Traj3DFileName))
        saveas(PlotFigure, ...
            fullfile(obj.SMF.Data.ResultsDir, [Traj3DFileName, '.png']))
        close(PlotFigure);
    elseif (obj.Verbose > 0)
        warning('smi.SPT.saveResults(): no plots produced, no trajectories!')
    end
end


end