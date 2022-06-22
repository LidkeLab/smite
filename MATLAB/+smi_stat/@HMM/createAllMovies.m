function [MovieParams] = createAllMovies(TRArray, SMFArray, ...
    MovieParams, SaveDir)
%createAllMovies creates dimer movies for all pairs in TRArray.

% INPUTS:
%   TRArray: TRArray with fields 'RawDataPath', 'DimerCandidateBool', and...
%            'StateSequence' populated.
%   SMFArray: Array of SMF structures corresponding to each trajectory pair
%             in TRArray.
%   MovieParams: Structure of movie parameters (see
%                smi_vis.GenerateMovies.prepDefaults()).
%                (Default = smi_vis.GenerateMovies.prepDefaults() with a
%                few updated settings)
%   SaveDir: Save directory for the movies. (Default = pwd())

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Set defaults.
if (~exist('SMFArray', 'var') || isempty(SMFArray))
    DefaultSMF = smi_core.SingleMoleculeFitting;
    SMFArray = repmat(DefaultSMF, size(TRArray, 1), 2);
end
if (~exist('MovieParams', 'var') || isempty(MovieParams))
    MovieParams = smi_vis.GenerateMovies.prepDefaults();
end
DefaultParams = smi_vis.GenerateMovies.prepDefaults();
DefaultParams.AutoCrop = true;
DefaultParams.IndicateDimer = true;
DefaultParams.IndicateDimerCandidate = true;
DefaultParams.TrajColors = [0, 1, 0; 1, 0, 1];
MovieParams = smi_helpers.padStruct(MovieParams, DefaultParams);
if (~exist('SaveDir', 'var') || isempty(SaveDir))
    SaveDir = pwd();
end
if ~exist(SaveDir, 'dir')
    mkdir(SaveDir);
end

% Loop through dimer pairs and prepare the movies.
for nn = 1:size(TRArray, 1)
    % Load the raw data for this dimer pair, if possible.
    FullFile1 = fullfile(SMFArray(nn, 1).Data.FileDir, ...
        SMFArray(nn, 1).Data.FileName{1});
    if isfile(FullFile1)
        LD = smi_core.LoadData;
        [~, RawData1] = LD.loadRawData(SMFArray(nn, 1), 1, ...
            SMFArray(nn, 1).Data.DataVariable);
    else
        RawData1 = ones(TRArray(nn, 1).YSize, ...
            TRArray(nn, 1).XSize);
    end
    FullFile2 = fullfile(SMFArray(nn, 2).Data.FileDir, ...
        SMFArray(nn, 2).Data.FileName{1});
    if isfile(FullFile2)
        LD = smi_core.LoadData;
        [~, RawData2] = LD.loadRawData(SMFArray(nn, 2), 1, ...
            SMFArray(nn, 2).Data.DataVariable);
    else
        RawData2 = ones(TRArray(nn, 2).YSize, ...
            TRArray(nn, 2).XSize);
    end

    % Load the registration transform, if available (we'll use this to
    % transform the raw data so that it aligns nicely for the movie).
    if isfile(SMFArray(nn, 2).Data.RegistrationFilePath)
        % We'll always assume channel 2 is the registered channel!
        load(SMFArray(nn, 2).Data.RegistrationFilePath, ...
            'RegistrationTransform')
        RawData2 = smi_core.ChannelRegistration.transformImages(...
            RegistrationTransform{2}, RawData2);
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