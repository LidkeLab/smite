% This script demonstrates the basic usage of the smi.SPT class.

%% Simulate and save some SPT data in the format expected for real data.
% Simulate some diffusing blobs.
rng(12)
SPTSim = smi_sim.SimSPT;
SPTSim.SimParams.FrameSize = [128, 128];
SPTSim.SimParams.ParticleDensity = 0.002; % particles / px^2
SPTSim.SimParams.D = 0.1; % px^2 / s
SPTSim.SimParams.KOffToOn = 0.9;
SPTSim.SimParams.KOnToOff = 0.05;
SPTSim.SimParams.Intensity = 1000;
SPTSim.createSimulation()
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.DataROI = [1, 1, SPTSim.SMD.YSize, SPTSim.SMD.XSize];
SMF.Fitting.PSFSigma = 1.3;
[~, sequence] = smi_sim.GaussBlobs.gaussBlobImage(SPTSim.SMD, SMF);

% Save the simulated data in a .mat file.
SaveDir = smi_helpers.mkSMITETmpDir('examples', 'SPT');
DataDir = fullfile(SaveDir, 'example_data', 'spt');
if ~isfolder(DataDir)
    mkdir(DataDir)
end
FileName = sprintf('Data_%s.mat', smi_helpers.genTimeString());
save(fullfile(DataDir, FileName), 'sequence')

%% Prepare a new SMF structure which points to the previously saved data.
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.FileDir = DataDir;
SMF.Data.FileName = FileName;
SMF.Fitting.PSFSigma = 1.3;
SMF.Thresholding.AutoThreshLogL = true;
% SMF.gui() % optionally, use SMF GUI to interactively set parameters

%% Prepare the SPT class object.
SPT = smi.SPT(SMF, false);
SPT.GenerateMovies = false;
SPT.GeneratePlots = false;

%% Perform the tracking.
SPT.performFullAnalysis()

%% Interactively prepare a movie.
% Reload our tracking results (this contains TR, SMD, and SMF).
[~, FileName] = fileparts(SMF.Data.FileName{1});
BaseName = [FileName, ...
    smi_helpers.arrayMUX({'_', ''}, isempty(SMF.Data.AnalysisID)), ...
    SMF.Data.AnalysisID];
ResultsFileName = [BaseName, '_Results.mat'];
ResultsFile = fullfile(DataDir, 'Results', ResultsFileName);
load(ResultsFile)

% Reload the raw data.
LD = smi_core.LoadData;
[~, RawData, SMF] = LD.loadRawData(SMF, 1, SMF.Data.DataVariable);

% Open up the movie maker.
MovieMaker = smi_vis.GenerateMovies;
MovieMaker.TR = TR;
MovieMaker.SMD = SMD;
MovieMaker.SMF = SMF;
MovieMaker.RawData = RawData;
MovieMaker.gui()



