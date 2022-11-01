
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
