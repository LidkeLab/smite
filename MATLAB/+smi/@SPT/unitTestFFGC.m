function [Success] = unitTestFFGC()
%unitTestFFGC tests the frame-to-frame and gap closing processes.
% This method was written with the intention of being a general smi.SPT
% testing script to test the internal SPT components (e.g., frame-to-frame
% connection, gap closing, and associated methods).  The idea behind this
% unit test is to test the tracking components without concern for the
% underlying localizations (i.e., we start from perfect localizations).


% Simulate Tracking Data (create a simulated SMD structure).
SPTSim = smi_sim.SimSPT;
SPTSim.SimParams.ParticleDensity = 0.01;
SPTSim.SimParams.Intensity = 500;
SPTSim.createSimulation()
SMD = SPTSim.SMD;
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.DataROI = [1, 1, SMD.YSize, SMD.XSize];
SMF.Fitting.PSFSigma = 1.3;
[~, RawData] = smi_sim.GaussBlobs.gaussBlobImage(SMD, SMF, 5);

% Perform the frame-to-frame connection process.
SMF = smi_core.SingleMoleculeFitting();
SMF.Tracking.D = SPTSim.SimParams.D * ones(size(SMD.FrameNum));
SMF.Tracking.K_off = SPTSim.SimParams.KOnToOff;
SMF.Tracking.K_on = SPTSim.SimParams.KOffToOn;
SMF.Tracking.MaxFrameGap = 20;
SMF.Tracking.MaxDistGC = 10;
SMF.Tracking.MaxDistFF = 2 * SMF.Tracking.MaxDistGC ...
    / sqrt(SMF.Tracking.MaxFrameGap);
SMD.ConnectID = zeros(numel(SMD.X), 1); 
RhoOnMean = mean(smi_core.SingleMoleculeData.computeDensity(SMD));
RhoOffMean = (SPTSim.SimParams.KOnToOff/SPTSim.SimParams.KOffToOn) * RhoOnMean;
SMF.Tracking.Rho_off = RhoOffMean;
for ff = min(SMD.FrameNum):(max(SMD.FrameNum)-1)
    [CM] = smi.SPT.createCostMatrixFF(SMD, SMF, RhoOffMean, ff, -1);
    [Link12] = smi.SPT.solveLAP(CM);
    [SMD] = smi.SPT.connectTrajFF(SMD, Link12, ff);
end

% Perform the gap-closing process.
[CM] =smi.SPT.createCostMatrixGC(SMD, SMF, RhoOffMean, -1, 1); 
[Link12] = smi.SPT.solveLAP(CM);
[SMD] = smi.SPT.connectTrajGC(SMD, Link12);
TR = smi_core.TrackingResults.convertSMDToTR(SMD);

% Make a movie.
MovieMaker = smi_vis.GenerateMovies;
MovieMaker.TR = TR;
MovieMaker.RawData = RawData;
MovieMaker.gui()

% Indicate success (this should be done in a better way, just setting it to
% 1 here isn't very useful).
Success = 1;


end