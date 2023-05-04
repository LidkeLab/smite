% This script demonstrates the basic usage of smi_stat.HMM.  We start by
% simulating some raw tracking data, reloading the raw data and tracking,
% and finally running the tracking results through the HMM analysis
% pipeline to search for dimers.  Simulated data and results will be saved
% in freshly generated directories under smite/MATLAB/examples/spt/

%% Simulate and save some SPT data in the format expected for real data.
% Simulate some diffusing blobs and save simulated raw data.
rng(12)
DataParams.NDatasets = 10;
SimParams.FrameSize = [64, 64];
SimParams.NFrames = 1000;
SimParams.ParticleDensity = 0.01; % particles / px^2, make 2X target since we split into 2 channels
SimParams.D = 0.1; % px^2 / frame
SimParams.KOffToOn = 0.95;
SimParams.KOnToOff = 0.05;
SimParams.KOnToBleach = 1e-3;
SimParams.Intensity = 1000;
SimParams.KDisconnect = 0.05;
SimParams.InteractionDistance = 0.5; % pixels
SimParams.InteractionProb = 0.5;
SaveDir = fullfile(tempdir, 'smite', 'examples', 'HMM');
if ~isfolder(SaveDir)
   mkdir(fullfile(tempdir, 'smite'));
   mkdir(fullfile(tempdir, 'smite', 'examples'));
   mkdir(fullfile(tempdir, 'smite', 'examples', 'HMM'));
end
DataDir = fullfile(SaveDir, ...
    'example_data', 'spt', smi_helpers.genTimeString());
if ~isfolder(DataDir)
    mkdir(DataDir)
end
[~, SimParams, DataParams] = smi_sim.SimSPT.makeExampleSim(...
    SimParams, DataParams, DataDir); % SimParams and DataParams padded w/ defaults

%% Prepare an SMF structure for each channel.
% NOTE: Channel 1 will not be transformed (it is the reference channel), so
%       we should ensure SMFChannel1.Data.RegistrationFilePath is empty!
PixelSize = 0.1; % must be set to true data values
FrameRate = 20;
SMFChannel1 = smi_core.SingleMoleculeFitting;
SMFChannel1.Data.AnalysisID = 'Channel1';
SMFChannel1.Data.FileDir = DataDir;
SMFChannel1.Data.DataROI = [1, 1, ...
    SimParams.FrameSize(1), SimParams.FrameSize(2)]; % left half of data [YStart, XStart, YEnd, XEnd]
SMFChannel1.Data.RegistrationFilePath = '';
SMFChannel1.Data.PixelSize = PixelSize;
SMFChannel1.Data.FrameRate = FrameRate;
SMFChannel1.Fitting.PSFSigma = SimParams.PSFSigma;
SMFChannel1.Tracking.MaxFrameGap = 20;
SMFChannel2 = copy(SMFChannel1);
SMFChannel2.Data.AnalysisID = 'Channel2';
SMFChannel2.Data.DataROI = [1, 1 + SimParams.FrameSize(2), ...
    SimParams.FrameSize(1), 2*SimParams.FrameSize(2)]; % right half of data
SMFChannel2.Data.RegistrationFilePath = '';

%% Prepare an SPT class object for each channel and then track.
SPTChannel1 = smi.SPT(SMFChannel1, false);
SPTChannel1.GenerateMovies = false;
SPTChannel1.GeneratePlots = false;
SPTChannel2 = smi.SPT(SMFChannel2, false);
SPTChannel2.GenerateMovies = false;
SPTChannel2.GeneratePlots = false;

%% Loop through our data files and track one at a time.
FileList = dir(fullfile(DataDir, 'Data_*.mat'));
FileList = FileList(~[FileList.isdir]); % exclude directories
NFiles = numel(FileList);
for nn = 1:NFiles
    % Update the file name in each of the SPT objects.
    SPTChannel1.SMF.Data.FileName = FileList(nn).name;
    SPTChannel2.SMF.Data.FileName = FileList(nn).name;
    
    % Perform the tracking.
    SPTChannel1.performFullAnalysis()
    SPTChannel2.performFullAnalysis()
end

%% Reload the tracked results and search for dimer candidates.
% NOTE: 'FileDir' is a directory containing both channel 1 and channel 2
%        results files.  findDimerCandidatesFromFiles() will search this
%        directory for files matching FilePatterns{1} and treat those files
%        as channel 1 results.  It will then do the same for
%        FilePatterns{2} (channel 2 results) and then attempt to match the
%        channel 1 and channel 2 files based on their file names (e.g.,
%        Data01_Channel1_Results.mat would be paired to
%        Data01_Channel2_Results.mat).
ResultsDir = fullfile(DataDir, 'Results');
FilePatterns = {'*Channel1_Results.mat'; '*Channel2_Results.mat'};
MaxDimerSeparation = 2; % pixels
MaxSeparation = 5; % pixels
[TRArray, ~, FileList] = smi_stat.HMM.findDimerCandidatesFromFiles(...
    ResultsDir, FilePatterns, MaxDimerSeparation, MaxSeparation);

% Add diffusion coefficients to the TRArray.
% NOTE: You can also do this on a per-trajectory basis if needed.
[TRArray.DiffusionCoefficient] = deal(SimParams.D);

%% Prepare the HMM class and run the analysis.
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.FrameRate = FrameRate; % fps, specific to the loaded TRs above
SMF.Data.PixelSize = PixelSize; % micrometers, specific to the loaded TRs above
HMM = smi_stat.HMM(TRArray, SMF);
HMM.MaxSeparation = MaxSeparation;
HMM.DimerSeparation = SimParams.InteractionDistance; % pixels, TRUE physical separation between dimers
HMM.SaveDir = fullfile(ResultsDir, 'HMM_Results');
HMM.GeneratePlots = [true; true]; % [basic plots; summary plots]
HMM.UnitFlag = false;
HMM.performFullAnalysis()
