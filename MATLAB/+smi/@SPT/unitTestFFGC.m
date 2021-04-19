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
SimParams.ParticleDensity = 0.01; % particles / pixel^2
SimParams.NFrames = 200;
SimParams.NSubFrames = 4;
SimParams.BoundaryCondition = 'Reflecting';
SimParams.Intensity = 1000; % photons / trajectory / frame
SimParams.D = 0.01; % pixel^2 / frame
SimParams.InteractionProb = 0;
SimParams.KBlinkOff = 0.3;
SimParams.KBlinkOn = 0.7;
SimParams.KBleach = 0;
SimParams.SigmaNoise = 0.2; % pixels
SimParams.FrameSize = 32; % pixels
SimParams.FrameRate = 1; % frames / second
SimParams.PixelSize = 0.1; % micrometers
SimParams.PSFSigma = 1.3; % pixels
SimParams.Bg = 5; % photons
[TD] = SMA_Sim.simulateTrajectories(SimParams);
TD.ConnectID = TD.TrajectoryID;
TD.NFrames = SimParams.NFrames;
[TRTruth] = smi_core.TrackingResults.convertSMDToTR(TD);
[RawData] = SMA_SPT.simRawDataFromTR(TRTruth, SimParams);

% Perform the frame-to-frame connection process.
SMF = smi_core.SingleMoleculeFitting();
SMF.Tracking.D = 0.01;
SMF.Tracking.K_off = 0.3;
SMF.Tracking.K_on = 0.7;
SMF.Tracking.MaxDistGC = 10;
SMF.Tracking.MaxDistFF = 2 * MaxDistGC / sqrt(MaxFrameGap);
SMF.Tracking.MaxFrameGap = 20;
TD.ConnectID = zeros(numel(TD.X), 1); 
RhoOnMean = mean(SMA_SPT.calcDensity(TD));
RhoOffMean = (SimParams.KBlinkOff/SimParams.KBlinkOn) * RhoOnMean;
SMF.Tracking.Rho_off = RhoOffMean;
for ff = min(TD.FrameNum):max(TD.FrameNum)
    [CM] = smi.SPT.createCostMatrixFF(TD, SMF, ff, -1);
    [Link12] = smi.SPT.solveLAP(CM);
    [TD] = smi.SPT.connectTrajFF(TD, Link12, ff);
end

% Perform the gap-closing process.
[CM] =smi.SPT.createCostMatrixGC(TD, SMF, -1, 1); 
[Link12] = smi.SPT.solveLAP(CM);
[TD] = smi.SPT.connectTrajGC(TD, Link12);
TR = smi_core.TrackingResults.convertSMDToTR(TD);

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