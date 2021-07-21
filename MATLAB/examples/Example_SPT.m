% This script demonstrates the basic usage of the smi.SPT class.

%% Prepare an SMF structure.
SMF = smi_core.SingleMoleculeFitting;
SMF.gui()

%% Prepare the SPT class object.
SPT = smi.SPT(SMF, false);
SPT.GenerateMovies = false;
SPT.GeneratePlots = false;

%% Perform the tracking.
SPT.performFullAnalysis()

%% Interactively prepare a movie.
% Reload our tracking results (this contains TR, SMD, and SMF).
load('C:\Users\David\Documents\GitHub\smite\MATLAB\examples\example_data\spt\Results\Data1__Results.mat')

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



