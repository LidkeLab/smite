
% This script demonstrates the basic usage of the smi_stat.ChangeDetection class.

%% simulate data
nObservations = 50;  % Scalar integer: length of data sequence
changePoints = round([0.2,0.6]*nObservations); %vector: indexs of change point locations
intensity = [200,100,60]; %vector: length=length(changePoints)+1: mean intensities for each subinterval
data=smi_stat.ChangeDetection.simulate(nObservations, changePoints, intensity);

% change detection
logBayesThreshold = 10;
icp = smi_stat.ChangeDetection(data,logBayesThreshold);
icp.plotIntensityEstimate();


%% alternative simulation methods
% given change points
nObservations = 50;  % Scalar integer: length of data sequence
changePoints = round([0.2,0.6]*nObservations); %vector: indexs of change point locations
intensity = [200,100,60]; %vector: length=length(changePoints)+1: mean intensities for each subinterval
logBayesThreshold = 10;
data=smi_stat.ChangeDetection.simulate(nObservations, changePoints, intensity);

[icp,f]=smi_stat.ChangeDetection.plotSimulatedEstimate(nObservations, changePoints, intensity, logBayesThreshold);

%% random change points
nObservations = 50;  % Scalar integer: length of data sequence
nChangePoints = 2; % Scalar integer: number of change points to simulate
meanintensity = 100; % scalar for mean intensity.  Intensities are uniformly distributed on [1, 2*meanIntensity]
logBayesThreshold = 10;
[icp,f]=smi_stat.ChangeDetection.plotRandSimulatedEstimate(nObservations, nChangePoints, meanintensity, logBayesThreshold);