% This script demonstrates the usage of smi_stat.DiffusionEstimator.

%% Single diffusing population, estimating D as the slope of the MSD.
% Simulate some trajectories.
SimParams = smi_sim.SimSPT.defineDefaultParams();
SimParams.ParticleDensity = 0.01; % particles / px^2
SimParams.NFrames = 1000;
SimParams.SubframeDensity = 1;
SimParams.FrameSize = [128, 128];
SimParams.BoundaryCondition = 'Periodic';
SimParams.KOffToOn = 1;
SimParams.KOnToOff = 0;
SimParams.D = 0.0123; % px^2 / frame
SPTSim = smi_sim.SimSPT(SimParams);
SPTSim.createSimulation();
TR = SPTSim.TR;

% Prepare the DiffusionEstimator class and estimate D.
% NOTE: For MSD fitting, the important parameters are 'DiffusionModel', 
%       'FitMethod', 'FrameLagRange', and 'NFitPoints'.  
SMF = smi_core.SingleMoleculeFitting; % used for PixelSize and FrameRate
DE = smi_stat.DiffusionEstimator(TR, SMF);
DE.FitTarget = 'MSD';
DE.DiffusionModel = 'Brownian'; % must be 'Brownian' for MSD fitting
DE.NComponents = 1; % must be 1 for MSD fitting
DE.FitMethod = 'WeightedLS'; % 'WeightedLS' or 'LS'
DE.FrameLagRange = [1, 50]; % range of frame lags computed in MSD
DE.NFitPoints = 5; % # of points in MSD used for fit.
DE.EstimateSEs = true; % estimate standard errors of fit D values
DE.FitIndividualTrajectories = true; % fit each trajectory MSD
DE.UnitFlag = false;
DE.Verbose = 1;
DE.estimateDiffusionConstant();

% Plot the MSD and fit results.
DE.plotEnsembleMSD(axes(figure()), DE.MSDEnsemble, DE.DiffusionStruct, ...
    DE.DiffusionModel, DE.UnitFlag);

%% Two diffusing populations, estimating D by fitting the CDF of jumps.
% Simulate some trajectories.
SimParams = smi_sim.SimSPT.defineDefaultParams();
SimParams.ParticleDensity = 0.01; % particles / px^2
SimParams.NFrames = 1000;
SimParams.SubframeDensity = 1;
SimParams.FrameSize = [128, 128];
SimParams.BoundaryCondition = 'Periodic';
SimParams.KOffToOn = 1;
SimParams.KOnToOff = 0;
SimParams.D = [0.0123, 0.00321, 0.00321]; % px^2 / frame
SPTSim = smi_sim.SimSPT(SimParams);
SPTSim.createSimulation();
TR = SPTSim.TR;

% Prepare the DiffusionEstimator class and estimate D.
% NOTE: For CDF fitting, the important parameters are 'DiffusionModel' and
%       'FrameLagRange'.
SMF = smi_core.SingleMoleculeFitting; % used for PixelSize and FrameRate
DE = smi_stat.DiffusionEstimator(TR, SMF);
DE.FitTarget = 'CDFOfJumps';
DE.DiffusionModel = 'Brownian';
DE.NComponents = 2;
DE.FrameLagRange = [1, 5]; % range of frame lags computed in MSD
DE.EstimateSEs = false; % estimate standard errors of fit D values
DE.FitIndividualTrajectories = false; % fit each trajectory MSD
DE.UnitFlag = false;
DE.Verbose = 1;
DE.estimateDiffusionConstant();

% Plot the CDF and fit results.
DE.plotEnsembleCDFOfJumps(axes(figure()), ...
    DE.MSDEnsemble, DE.DiffusionStruct, DE.UnitFlag);

%% Two diffusing populations, estimating D with an MLE for the jumps.
% (If estimating standard errors, this is faster than fitting the CDF, and
% the results are similar).

% Prepare the DiffusionEstimator class and estimate D.
% NOTE: For LikelihoodOfJumps maximization, the important parameters are
%       'DiffusionModel' and 'FrameLagRange'.
SMF = smi_core.SingleMoleculeFitting; % used for PixelSize and FrameRate
DE = smi_stat.DiffusionEstimator(TR, SMF);
DE.FitTarget = 'LikelihoodOfJumps';
DE.DiffusionModel = 'Brownian';
DE.NComponents = 2;
DE.FrameLagRange = [1, 5]; % range of frame lags computed in MSD
DE.EstimateSEs = true; % estimate standard errors of fit D values
DE.FitIndividualTrajectories = false; % fit each trajectory MSD
DE.UnitFlag = false;
DE.Verbose = 1;
DE.estimateDiffusionConstant();
