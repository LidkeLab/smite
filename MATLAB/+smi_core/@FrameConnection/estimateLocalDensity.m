function [InitialDensity] = ...
    estimateLocalDensity(ClusterData, NNearestClusters, ...
    KOn, KOff, KBleach, PMiss)
%estimateLocalDensity estimates the local density for clustered locs.
% This method will make an estimate of the local density of clusters for
% each cluster in ClusterData.
%
% INPUTS:
%   ClusterData: Cell array of cluster data (see organizeClusterData())
%   NNearestClusters: Number of nearest clusters that we'll search for to
%                     estimate the density.
%   KOn: Rate parameter for blinking on. (1 / frame)
%   KOff: Rate parameter for blinking off. (1 / frame)
%   KBleach: Rate parameter for photobleaching. (1 / frame)
%   PMiss: Probability of missing a localization of a visible emitter. 
%
% OUTPUTS:
%   InitialDensity: Estimate of underlying emitter density at the start of
%                   the experiment for each cluster. (emitters / pixel^2)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Compute the center of each cluster.  If there is only one cluster we'll
% stop early before doing the knnsearch() below.
DutyCycle = KOn / (KOn+KOff+KBleach);
NClusters = numel(ClusterData);
MaxFrame = double(max(cellfun(@(X) max(X(:, 5) .* X(:, 6)), ClusterData)));
if (NClusters == 1)
    % If we only have 1 cluster, we'll just assume it's spatial extent is
    % the area.
    if (numel(ClusterData{1}(:, 1)) == 1)
        InitialDensity = 1;
    else
        Area = (max(ClusterData{1}(:, 1))-min(ClusterData{1}(:, 1))) ...
            * (max(ClusterData{1}(:, 2))-min(ClusterData{1}(:, 2)));
        InitialDensity = (1/Area) * ((KBleach/KOff) / (1-PMiss)) ...
            / (1-exp(-KBleach*DutyCycle*(MaxFrame-1)));
    end
    return
end
ClusterCenters = NaN(NClusters, 2);
for nn = 1:NClusters
    % Isolate some arrays to improve readability.
    X = ClusterData{nn}(:, 1);
    Y = ClusterData{nn}(:, 2);
    X_SE = ClusterData{nn}(:, 3);
    Y_SE = ClusterData{nn}(:, 4);
    
    % Compute the cluster center to be the MLE of the position if all
    % cluster localizations are combined into one.
    ClusterCenters(nn, :) = [sum(X./X_SE.^2) / sum(1./X_SE.^2), ...
        sum(Y./Y_SE.^2) / sum(1./Y_SE.^2)];
end

% For each cluster, estimate the cluster density by assuming each of the
% NNearestClusters are individual emitters.
NNearestClusters = min(NNearestClusters, NClusters-1);
[~, NNClusterDistances] = knnsearch(ClusterCenters, ClusterCenters, ...
    'k', NNearestClusters + 1);
LocalClusterDensity = (NNearestClusters+1) ...
    ./ (pi*NNClusterDistances(:, NNearestClusters+1).^2);

% Based on the cluster density, make an estimate of the initial underlying
% emitter density.
Lambda1 = KBleach * DutyCycle;
Lambda2 = (KOn+KOff+KBleach) - Lambda1;
InitialDensity = LocalClusterDensity ...
    * (1/DutyCycle) * (1/KOff) * (1/(1-PMiss)) ...
    ./ ((1/Lambda1)*(1-exp(-Lambda1*(MaxFrame-1))) ...
    - (1/Lambda2)*(1-exp(-Lambda2*(MaxFrame-1))));


end