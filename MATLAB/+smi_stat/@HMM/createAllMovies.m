function [MovieParams] = createAllMovies(TRArray, ...
    MovieParams, SaveDir, RawDataBaseDir)
%createAllMovies creates dimer movies for all pairs in TRArray.

% INPUTS:
%   TRArray: TRArray with fields 'RawDataPath', 'DimerCandidateBool', and...
%            'StateSequence' populated.
%   MovieParams: Structure of movie parameters (see
%                smi_vis.GenerateMovies.prepDefaults()).
%                (Default = smi_vis.GenerateMovies.prepDefaults() with a
%                few updated settings)
%   SaveDir: Save directory for the movies. (Default = pwd())
%   RawDataBaseDir: Raw data directory.  This can be useful to use a 
%                   different path than the one saved in 'RawDataPath'.
%                   (Default chosen from data in TRArray)

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Set defaults.
if (~exist('MovieParams', 'var') || isempty(MovieParams))
    MovieParams = smi_vis.GenerateMovies.prepDefaults();
end
DefaultParams = smi_vis.GenerateMovies.prepDefaults();
DefaultParams.AutoCrop = true;
DefaultParams.IndicateDimer = true;
DefaultParams.IndicateDimerCandidate = true;
MovieParams = smi_helpers.padStruct(MovieParams, DefaultParams);
if (~exist('SaveDir', 'var') || isempty(SaveDir))
    SaveDir = pwd();
end
if ~exist(SaveDir, 'dir')
    mkdir(SaveDir);
end
if (~exist('RawDataBaseDir', 'var') || isempty(RawDataBaseDir))
    RawDataBaseDir = fileparts(TRArray(1).RawDataPath);
end


% Loop through dimer pairs and prepare the movies.
SMF = smi_core.SingleMoleculeFitting;
for nn = 1:size(TRArray, 1)
    % Load the raw data for this dimer pair, if possible.
    if isfile(TRArray(nn, 1).RawDataPath)
        [~, File, Extension] = ...
            fileparts(TRArray(nn, 1).RawDataPath);
        SMF.Data.FileDir = RawDataBaseDir;
        SMF.Data.FileName = [File, Extension];
        SMF.Data.FileType = Extension;
        LD = smi_core.LoadData;
        [~, RawData1] = LD.loadRawData(SMF, 1, ...
            SMF.Data.DataVariable);
    else
        RawData1 = ones(TRArray(nn, 1).YSize, ...
            TRArray(nn, 1).XSize);
    end
    if isfile(TRArray(nn, 2).RawDataPath)
        [~, File, Extension] = ...
            fileparts(TRArray(nn, 2).RawDataPath);
        SMF.Data.FileDir = RawDataBaseDir;
        SMF.Data.FileName = [File, Extension];
        SMF.Data.FileType = Extension;
        LD = smi_core.LoadData;
        [~, RawData2] = LD.loadRawData(SMF, 1, ...
            SMF.Data.DataVariable);
    else
        RawData2 = ones(TRArray(nn, 1).YSize, ...
            TRArray(nn, 1).XSize);
    end

    % Rescale the raw data and create the false color overlay.
    RawData1 = smi_vis.contrastStretch(RawData1, [], ...
        MovieParams.PercentileCeiling, ...
        MovieParams.PercentileFloor, ...
        MovieParams.MinScaleIntensity);
    RawData2 = smi_vis.contrastStretch(RawData2, [], ...
        MovieParams.PercentileCeiling, ...
        MovieParams.PercentileFloor, ...
        MovieParams.MinScaleIntensity);
    RawDataRGB = ones([size(RawData1, 1:2), 3, size(RawData1, 3)]);
    for cc = 1:3
        RawDataRGB(:, :, cc, :) = ...
            RawData1*MovieParams.RawDataColors(1, cc) ...
            + RawData2*MovieParams.RawDataColors(2, cc);
    end

    % Prepare and save the movie.
    MovieMaker = smi_vis.GenerateMovies(MovieParams);
    MovieMaker.RawData = RawDataRGB;
    MovieMaker.TR = TRArray(nn, :);
    MovieSavePath = fullfile(SaveDir, ...
        sprintf('DimerPair%i_movie.mp4', nn));
    MovieMaker.saveMovie(MovieSavePath)
end


end