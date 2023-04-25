% Test calling smi_stat.ChangeDetection routines.
function success = unitTest()

success = 0;

% This script demonstrates the basic usage of the smi_stat.ChangeDetection class.

%% simulate data
NObservations = 50;  % Scalar integer: length of data sequence
ChangePoints = round([0.2,0.6]*NObservations); %vector: indexs of change point locations
Intensity = [200,100,60]; %vector: length=length(ChangePoints)+1: mean intensities for each subinterval
Data=smi_stat.ChangeDetection.simulate(NObservations, ChangePoints, Intensity);

% change detection
LogBayesThreshold = 10;
Icp = smi_stat.ChangeDetection(Data,LogBayesThreshold);
Icp.plotIntensityEstimate();
IntensityModel = Icp.IntensityModel; % intensity of the model sequence (no noise)

%% alternative simulation methods
% given change points
NObservations = 50;  % Scalar integer: length of data sequence
ChangePoints = round([0.2,0.6]*NObservations); %vector: indexs of change point locations
Intensity = [200,100,60]; %vector: length=length(ChangePoints)+1: mean intensities for each subinterval
LogBayesThreshold = 10;
Data=smi_stat.ChangeDetection.simulate(NObservations, ChangePoints, Intensity);

[Icp,F]=smi_stat.ChangeDetection.plotSimulatedEstimate(NObservations, ChangePoints, Intensity, LogBayesThreshold);

%% random change points
NObservations = 50;  % Scalar integer: length of data sequence
NChangePoints = 2; % Scalar integer: number of change points to simulate
Meanintensity = 100; % scalar for mean intensity.  Intensities are uniformly distributed on [1, 2*meanIntensity]
LogBayesThreshold = 10;
[Icp,F]=smi_stat.ChangeDetection.plotRandSimulatedEstimate(NObservations, NChangePoints, Meanintensity, LogBayesThreshold);

success = 1;

end
