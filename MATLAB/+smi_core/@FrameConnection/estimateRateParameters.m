function [KOn, KOff, KBleach, PMiss, NEmitters] = ...
    estimateRateParameters(ClusterData, Verbose)
%estimateRateParameters estimates rates from clustered localizations.
% This method will make an estimate of the blinking kinetics (KOn, KOff,
% and KBleach) based on the (pre-clustered) localizations in ClusterData.
% These are computed by assuming that each cluster corresponds to a single
% emitter.
%
% INPUTS:
%   ClusterData: Cell array of cluster data (see organizeClusterData())
%
% OUTPUTS:
%   KOn: Rate parameter for blinking on. (1 / frame)
%   KOff: Rate parameter for blinking off. (1 / frame)
%   KBleach: Rate parameter for photobleaching. (1 / frame)
%   PMiss: Probability of missing a localization.
%   NEmitters: Total number of emitters present at the start of the
%              experiment.
%   Verbose: Flag to indicate verbosity of outputs. (Default = 1)
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('Verbose', 'var') || isempty(Verbose))
    Verbose = 1;
end    

% Compute some quantities needed from each cluster.
NClusters = numel(ClusterData);
TotalDurations = NaN(NClusters, 1);
NObservations = NaN(NClusters, 1);
for nn = 1:NClusters
    % Compute the total duration of the cluster.
    CurrentFrames = ClusterData{nn}(:, 5);
    TotalDurations(nn) = max(CurrentFrames) - min(CurrentFrames) + 1;
    
    % Compute the total number of observed localizations (clusters might
    % have multiple localizations per frame due to generous pre-clustering,
    % so we should make sure not to overcount those).
    NObservations(nn) = numel(unique(CurrentFrames));
end

% Estimate KOff+KBleach and PMiss, assuming each cluster was a single 
% blinking event of a single emitter.
KOffPKBleach = -log(1 - 1/mean(TotalDurations));
PMiss = mean(abs((NObservations-TotalDurations)) ./ TotalDurations);

% Compute some parameters from the sum of localizations present over time.
AllData = cell2mat(ClusterData);
[NLoc, Frames] = groupcounts(double(AllData(:, 5)));
NLocSum = cumsum(NLoc);
CostFunction = @(X) mean((NLocSum - X(1)*(1-exp(-X(2)*(Frames-1)))).^2);
FitOptions = optimset('Display', ...
    smi_helpers.stringMUX({'none', 'final'}, (Verbose>1)));
LocSumParams = fmincon(CostFunction, ...
    [NLocSum(end), 1/Frames(end)], [], [], [], [], [0, 0], [], [], ...
    FitOptions);
KOff = NClusters / NLocSum(end);
KBleach = KOffPKBleach - KOff;
NEmitters = LocSumParams(1) * KBleach / (1-PMiss);
KOn = LocSumParams(2) * KOffPKBleach / (KBleach-LocSumParams(2));


end