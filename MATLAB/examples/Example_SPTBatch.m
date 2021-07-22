% This script demonstrates an example workflow for SPT batch processing
% including the use of channel registration.

%% Prepare the channel registration transform.
FiducialDir = 'C:\Users\David\Documents\GitHub\smite\MATLAB\examples\example_data\spt\multiple_files\fiducial';
FiducialFile = {'Fiducial.mat'};
SMFFiducial = smi_core.SingleMoleculeFitting;
SMFFiducial.Fitting.FitType = 'XYNBS';
Verbose = 1;
ChannelReg = smi_core.ChannelRegistration(...
    FiducialDir, FiducialFile, SMFFiducial, Verbose);
ChannelReg.SplitFormat = [1, 2]; % split fiducial by columns at the center
ChannelReg.AutoscaleFiducials = true;
ChannelReg.ManualCull = true;
ChannelReg.TransformationType = 'lwm';
ChannelReg.TransformationBasis = 'coordinates';
ChannelReg.NNeighborPoints = 12;
ChannelReg.findTransform();
TransformPath = ChannelReg.exportTransform();

%% Prepare an SMF structure for each channel.
% NOTE: Channel 1 will not be transformed (it is the reference channel), so
%       we should ensure SMFChannel1.Data.RegistrationFilePath is empty!
FileDir = 'C:\Users\David\Documents\GitHub\smite\MATLAB\examples\example_data\spt\multiple_files';
PixelSize = 0.1;
FrameRate = 20;
SMFChannel1 = smi_core.SingleMoleculeFitting;
SMFChannel1.Data.AnalysisID = 'Channel1';
SMFChannel1.Data.FileDir = FileDir;
SMFChannel1.Data.DataROI = [1, 1, 64, 64]; % left half of data
SMFChannel1.Data.RegistrationFilePath = '';
SMFChannel1.Data.PixelSize = PixelSize;
SMFChannel1.Data.FrameRate = FrameRate;
SMFChannel1.Fitting.PSFSigma = 1.3;
SMFChannel1.Tracking.MaxFrameGap = 20;
SMFChannel2 = smi_core.SingleMoleculeFitting;
SMFChannel2.Data.AnalysisID = 'Channel2';
SMFChannel2.Data.FileDir = FileDir;
SMFChannel2.Data.DataROI = [1, 65, 64, 128]; % right half of data
SMFChannel2.Data.RegistrationFilePath = TransformPath{1};
SMFChannel2.Fitting.PSFSigma = 1.3;
SMFChannel2.Tracking.MaxFrameGap = 20;
SMFChannel2.Data.PixelSize = PixelSize;
SMFChannel2.Data.FrameRate = FrameRate;

%% Prepare an SPT class object for each channel and then track.
SPTChannel1 = smi.SPT(SMFChannel1, false);
SPTChannel1.GenerateMovies = false;
SPTChannel1.GeneratePlots = true;
SPTChannel2 = smi.SPT(SMFChannel2, false);
SPTChannel2.GenerateMovies = false;
SPTChannel2.GeneratePlots = true;

%% Loop through our data files and track one at a time.
FileList = dir(fullfile(FileDir, 'Data_*.mat'));
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

%% Interactively prepare a movie for each channel.
% Define the index of 'FileList' which we'd like to see a movie for.
FileIndex = 7;

% Reload our tracking results (these contain TR, SMD, and SMF).
[~, Channel1FileName] = fileparts(FileList(FileIndex).name);
Channel1ResultsFiles = dir(fullfile(SMFChannel1.Data.ResultsDir, ...
    sprintf('*%s_%s_Results.mat', ...
    Channel1FileName, SMFChannel1.Data.AnalysisID)));
Channel1Results = ...
    load(fullfile(SMFChannel1.Data.ResultsDir, Channel1ResultsFiles(1).name));
[~, Channel2FileName] = fileparts(FileList(FileIndex).name);
Channel2ResultsFiles = dir(fullfile(SMFChannel2.Data.ResultsDir, ...
    sprintf('*%s_%s_Results.mat', ...
    Channel2FileName, SMFChannel2.Data.AnalysisID)));
Channel2Results = ...
    load(fullfile(SMFChannel2.Data.ResultsDir, Channel2ResultsFiles(1).name));

% Reload the raw data.
SMFChannel1.Data.FileName = FileList(FileIndex).name;
SMFChannel2.Data.FileName = FileList(FileIndex).name;
LD = smi_core.LoadData;
[~, RawDataChannel1] = ...
    LD.loadRawData(SMFChannel1, 1, SMFChannel1.Data.DataVariable);
[~, RawDataChannel2] = ...
    LD.loadRawData(SMFChannel2, 1, SMFChannel2.Data.DataVariable);

% If needed, transform the raw data (to make sure the trajectories overlay
% nicely on the raw data!).
if ~isempty(SMFChannel2.Data.RegistrationFilePath)
    load(SMFChannel2.Data.RegistrationFilePath, 'RegistrationTransform')
end
RawDataChannel2 = smi_core.ChannelRegistration.transformImages(...
    RegistrationTransform, RawDataChannel2);

% Prepare the movies using the GUI.
MovieMaker1 = smi_vis.GenerateMovies;
MovieMaker1.TR = Channel1Results.TR;
MovieMaker1.SMD = Channel1Results.SMD;
MovieMaker1.SMF = Channel1Results.SMF;
MovieMaker1.RawData = RawDataChannel1;
MovieMaker1.gui()
MovieMaker2 = smi_vis.GenerateMovies;
MovieMaker2.TR = Channel2Results.TR;
MovieMaker2.SMD = Channel2Results.SMD;
MovieMaker2.SMF = Channel2Results.SMF;
MovieMaker2.RawData = RawDataChannel2;
MovieMaker2.gui()



