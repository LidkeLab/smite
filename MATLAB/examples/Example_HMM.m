% Simulate some trajectories.
% NOTE: For now, this requires sma-core-alpha to use
%       SMA_Sim.simulateTrajectories(), which hasn't been ported to smite.
SimParams.AllowOligomers = 0; % allows higher order oligomers
SimParams.ParticleDensity = 0.005; % particles / px^2
SimParams.NSubFrames = 1000;
SimParams.NFrames = 1000;
SimParams.FrameSize = 128; % pixels (square only)
SimParams.PixelSize = 0.1;
SimParams.FrameRate = 20;
SimParams.InteractionDistance = 0.5;
SimParams.InteractionProb = 0.5;
SimParams.BoundaryCondition = 'Periodic';
SimParams.KBlinkOn = 1;
SimParams.KBlinkOff = 0;
SimParams.D = 0.0123 ...
    / (SimParams.FrameRate * SimParams.PixelSize^2); % px^2 / subframe
SimParams.Intensity = 300;
[SMD, ~, ~, SimParams] = SMA_Sim.simulateTrajectories(SimParams);
SMD.ConnectID = SMD.TrajectoryID; % fix for SMA_Sim.simulateTrajectories()
TR = smi_core.TrackingResults.convertSMDToTR(SMD);
TR1 = TR(1:floor(numel(TR)/2));
TR2 = TR(floor(numel(TR)/2)+1:end);

%% 
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.FrameRate = SimParams.FrameRate;
SMF.Data.PixelSize = SimParams.PixelSize;
MaxDimerSeparation = 2;
MaxSeparation = 5;
TRArray = smi_stat.HMM.findDimerCandidates(TR1, TR2, ...
    MaxDimerSeparation, MaxSeparation);
HMM = smi_stat.HMM(TRArray, SMF);
HMM.DiffusionCoefficient = SimParams.D;
HMM.DimerSeparation = SimParams.InteractionDistance;
HMM.MaxSeparation = MaxSeparation;
HMM.performFullAnalysis()




