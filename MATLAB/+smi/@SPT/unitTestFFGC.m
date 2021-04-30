function [Success] = unitTestFFGC()
%unitTestFFGC tests the frame-to-frame and gap closing processes.
% This method was written with the intention of being a general smi.SPT
% testing script to test the internal SPT components (e.g., frame-to-frame
% connection, gap closing, and associated methods).  The idea behind this
% unit test is to test the tracking components without concern for the
% underlying localizations (i.e., we start from perfect localizations).
% For now, we're still relying on some older sma-core-alpha methods to make
% use of this unit test.


% Simulate Tracking Data (create a simulated TD structure).
SimParams.ParticleDensity = 0.001; % particles / pixel^2
SimParams.NFrames = 500;
SimParams.NSubFrames = 1;
SimParams.BoundaryCondition = 'Free';
SimParams.Intensity = 1000; % photons / trajectory / frame
SimParams.D = 0.01; % pixel^2 / frame
SimParams.InteractionProb = 0;
SimParams.KBlinkOff = 0.2;
SimParams.KBlinkOn = 0.7;
SimParams.KBleach = 0;
SimParams.FrameSize = 64; % pixels
SimParams.FrameRate = 1; % frames / second
SimParams.PixelSize = 0.1; % micrometers
SimParams.PSFSigma = 1.3; % pixels
SimParams.SigmaNoise = ...
    SimParams.PSFSigma / sqrt(SimParams.Intensity); % pixels
SimParams.Bg = 5; % photons
[TD] = SMA_Sim.simulateTrajectories(SimParams);
TD.ConnectID = TD.TrajectoryID;
TD.NFrames = SimParams.NFrames;
[TRTruth] = smi_core.TrackingResults.convertSMDToTR(TD);
[TRTruth.PixelSize] = deal(SimParams.PixelSize);
[TRTruth.FrameRate] = deal(SimParams.FrameRate);
[RawData] = SMA_SPT.simRawDataFromTR(TRTruth, SimParams);

% Perform the frame-to-frame connection process.
SMF = smi_core.SingleMoleculeFitting();
SMF.Tracking.D = SimParams.D;
SMF.Tracking.K_off = SimParams.KBlinkOff;
SMF.Tracking.K_on = SimParams.KBlinkOn;
SMF.Tracking.MaxFrameGap = 20;
SMF.Tracking.MaxDistGC = 10;
SMF.Tracking.MaxDistFF = 2 * SMF.Tracking.MaxDistGC ...
    / sqrt(SMF.Tracking.MaxFrameGap);
TD.ConnectID = zeros(numel(TD.X), 1); 
RhoOnMean = mean(SMA_SPT.calcDensity(TD));
RhoOffMean = (SimParams.KBlinkOff/SimParams.KBlinkOn) * RhoOnMean;
SMF.Tracking.Rho_off = RhoOffMean;
for ff = min(TD.FrameNum):(max(TD.FrameNum)-1)
    [CM] = smi.SPT.createCostMatrixFF(TD, SMF, [], ff, -1);
    [Link12] = smi.SPT.solveLAP(CM);
    [TD] = smi.SPT.connectTrajFF(TD, Link12, ff);
end

% Perform the gap-closing process.
[CM] =smi.SPT.createCostMatrixGC(TD, SMF, [], -1, 1); 
[Link12] = smi.SPT.solveLAP(CM);
[TD] = smi.SPT.connectTrajGC(TD, Link12);
TR = smi_core.TrackingResults.convertSMDToTR(TD);
[TR.PixelSize] = deal(SimParams.PixelSize);
[TR.FrameRate] = deal(SimParams.FrameRate);

% Make some plots.
PlotFigure = figure();
SMA_SPT.plot2D(PlotFigure, TR)
PlotFigure = figure();
SMA_SPT.plot3D(PlotFigure, TR)
PlotFigure = figure();
DisplayParams.MaxTrajDisplayLength = 10; % frames
DisplayParams.AutoPlay = 0;
SMA_SPT.movieTraj(PlotFigure, TR, RawData, [], [], DisplayParams);

% Indicate success (this should be done in a better way, just setting it to
% 1 here isn't very useful).
Success = 1;


end