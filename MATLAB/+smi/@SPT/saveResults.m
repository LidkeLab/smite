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

% Create the filenames for the saved results.
% NOTE: For now, this is only setup for a single file (hence the {1}).
[~, FileName] = fileparts(obj.SMF.Data.FileName{1});
SMDFileName = sprintf('%s_%s_%s.mat', ...
    FileName, obj.SMF.Data.AnalysisID, 'SMD');
TRFileName = sprintf('%s_%s_%s.mat', ...
    FileName, obj.SMF.Data.AnalysisID, 'TR');
SMFFileName = sprintf('%s_%s_%s.mat', ...
    FileName, obj.SMF.Data.AnalysisID, 'SMF');

% Move the data structures of interest into the workspace with appropriate
% names for saving.
SMD = obj.SMDPreThresh;
TR = obj.TR;
SMF = obj.SMF;

% Save the data structures as .mat files.
save(fullfile(obj.SMF.Data.ResultsDir, SMDFileName), 'SMD');
save(fullfile(obj.SMF.Data.ResultsDir, TRFileName), 'TR');
save(fullfile(obj.SMF.Data.ResultsDir, SMFFileName), 'SMF');

% Create a movie of the tracks and save the resulting movie.
if obj.GenerateMovies
    % Load the raw data.
    MovieFileName = sprintf('%s_%s_movie.mp4', ...
        FileName, obj.SMF.Data.AnalysisID);
    LD = smi_core.LoadData;
    [~, RawData, obj.SMF] = ...
        LD.loadRawData(obj.SMF, 1, obj.SMF.Data.DataVariable);
    
    % Generate the movie.
    % NOTE: For now, we're using the OLD version of the movie code in
    %       SMA_SPT (we'll need to change this once a new version is
    %       written for smite).
    DisplayParams.UnitFlag = obj.UnitFlag;
    DisplayParams.LiteMode = 1;
    PlotFigure = figure();
    SMA_SPT.movieTraj(PlotFigure, obj.TR, RawData, [], ...
        fullfile(obj.SMF.Data.ResultsDir, MovieFileName), DisplayParams);
    close(PlotFigure)
end

% Create the 2D and 3D tracking results and save these figures.
if obj.GeneratePlots
    PlotFigure = figure();
    SMA_SPT.plot2D(PlotFigure, obj.TR, obj.UnitFlag, ...
        1, obj.SMF.Data.ResultsDir);
    clf(PlotFigure);
    SMA_SPT.plot3D(PlotFigure, obj.TR, obj.UnitFlag, ...
        1, obj.SMF.Data.ResultsDir);
    close(PlotFigure);
end


end