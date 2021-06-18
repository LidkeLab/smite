function [ClusterData] = organizeClusterData(SMD)
%organizeClusterData organizes data according to SMD.ConnectID.
% This method organizes relevant data (X, Y, ...) for each cluster into a
% cell of cell array for later use.
%
% INPUTS:
%   SMD: Single Molecule Data structure with a populated SMD.ConnectID.
%
% OUTPUTS:
%   ClusterData: Cell array of cluster data, with the index corresponding
%                to the sorted unique cluster IDs (e.g., if 
%                unique(ConnectID) = [1, 2, 5], 
%                ClusterData{1}<->ConnectID 1, 
%                ClusterData{2}<->ConnectID 2,
%                ClusterData{3}<->ConnectID 5

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Isolate/organize some SMD arrays.
ConnectID = SMD.ConnectID;
CombinedData = [SMD.X, SMD.Y, SMD.X_SE, SMD.Y_SE, ...
    single(SMD.FrameNum), single(SMD.DatasetNum), single(ConnectID)];
[~, SortIndices] = sort(ConnectID);
CombinedData = CombinedData(SortIndices, :);

% Loop over the unique connect IDs and prepare the output ClusterData.
NLocPerID = groupcounts(ConnectID);
NLocCumulative = [0; cumsum(NLocPerID)];
UniqueIDs = unique(ConnectID);
NUnique = numel(UniqueIDs);
ClusterData = cell(NUnique, 1);
for ii = 1:NUnique
    IndexArray = (1:NLocPerID(ii)) + NLocCumulative(ii);
    ClusterData{ii} = [CombinedData(IndexArray, :), ...
        SortIndices(IndexArray.')];
end


end