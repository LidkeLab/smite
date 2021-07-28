% This script demonstrates the usage of smi_stat.DiffusionEstimator.

%% Single diffusing population, estimating D as the slope of the MSD.
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
SimParams.BoundaryCondition = 'Periodic';
SimParams.KBlinkOn = 1;
SimParams.KBlinkOff = 0;
SimParams.D = 0.0123 ...
    / (SimParams.FrameRate * SimParams.PixelSize^2); % px^2 / subframe
SimParams.Intensity = 300;
[SMD, ~, ~, SimParams] = SMA_Sim.simulateTrajectories(SimParams);
SMD.ConnectID = SMD.TrajectoryID; % fix for SMA_Sim.simulateTrajectories()
TR = smi_core.TrackingResults.convertSMDToTR(SMD);

% Prepare the DiffusionEstimator class and estimate D.
% NOTE: For MSD fitting, the important parameters are 'DiffusionModel', 
%       'FitMethod', 'FrameLagRange', and 'NFitPoints'.  
SMF = smi_core.SingleMoleculeFitting; % used for PixelSize and FrameRate
SMF.Data.PixelSize = SimParams.PixelSize;
SMF.Data.FrameRate = SimParams.FrameRate;
DE = smi_stat.DiffusionEstimator(TR, SMF);
DE.FitTarget = 'MSD';
DE.DiffusionModel = 'Brownian1C'; % must be 'Brownian1C' for MSD fitting
DE.FitMethod = 'WeightedLS'; % 'WeightedLS' or 'LS'
DE.FrameLagRange = [1, 50]; % range of frame lags computed in MSD
DE.NFitPoints = 5; % # of points in MSD used for fit.
DE.EstimateSEs = true; % estimate standard errors of fit D values
DE.FitIndividualTrajectories = true; % fit each trajectory MSD
DE.UnitFlag = true;
DE.Verbose = 1;
DE.estimateDiffusionConstant();

% Plot the MSD and fit results.
DE.plotEnsembleMSD(axes(figure()), DE.MSDEnsemble, DE.DiffusionStruct, ...
    DE.DiffusionModel, DE.UnitFlag);

%% Two diffusing populations, estimating D by fitting the CDF of jumps.
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
SimParams.BoundaryCondition = 'Periodic';
SimParams.KBlinkOn = 1;
SimParams.KBlinkOff = 0;
SimParams.D = [0.0123, 0.00321, 0.00321] ...
    / (SimParams.FrameRate * SimParams.PixelSize^2); % px^2 / subframe
SimParams.Intensity = 300;
[SMD, ~, ~, SimParams] = SMA_Sim.simulateTrajectories(SimParams);
SMD.ConnectID = SMD.TrajectoryID; % fix for SMA_Sim.simulateTrajectories()
TR = smi_core.TrackingResults.convertSMDToTR(SMD);

% Prepare the DiffusionEstimator class and estimate D.
% NOTE: For CDF fitting, the important parameters are 'DiffusionModel' and
%       'FrameLagRange'.
SMF = smi_core.SingleMoleculeFitting; % used for PixelSize and FrameRate
SMF.Data.PixelSize = SimParams.PixelSize;
SMF.Data.FrameRate = SimParams.FrameRate;
DE = smi_stat.DiffusionEstimator(TR, SMF);
DE.FitTarget = 'CDFOfJumps';
DE.DiffusionModel = 'Brownian2C';
DE.FrameLagRange = [2, 2]; % range of frame lags computed in MSD
DE.EstimateSEs = false; % estimate standard errors of fit D values
DE.FitIndividualTrajectories = false; % fit each trajectory MSD
DE.UnitFlag = true;
DE.Verbose = 1;
DE.estimateDiffusionConstant();

% Plot the CDF and fit results.
DE.plotEnsembleCDFOfJumps(axes(figure()), ...
    DE.MSDEnsemble, DE.DiffusionStruct, DE.DiffusionModel, DE.UnitFlag);

%% Two diffusing populations, estimating D with an MLE for the jumps.
% (If estimating standard errors, this is faster than fitting the CDF, and
% the results are similar).

% Prepare the DiffusionEstimator class and estimate D.
% NOTE: For LikelihoodOfJumps maximization, the important parameters are
%       'DiffusionModel' and 'FrameLagRange'.
SMF = smi_core.SingleMoleculeFitting; % used for PixelSize and FrameRate
SMF.Data.PixelSize = SimParams.PixelSize;
SMF.Data.FrameRate = SimParams.FrameRate;
DE = smi_stat.DiffusionEstimator(TR, SMF);
DE.FitTarget = 'LikelihoodOfJumps';
DE.DiffusionModel = 'Brownian2C';
DE.FrameLagRange = [2, 2]; % range of frame lags computed in MSD
DE.EstimateSEs = true; % estimate standard errors of fit D values
DE.FitIndividualTrajectories = false; % fit each trajectory MSD
DE.UnitFlag = true;
DE.Verbose = 1;
DE.estimateDiffusionConstant();