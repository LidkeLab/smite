function [ConnectID, MaxConnectID] = ...
    linkClusters(ConnectID, MaxConnectID, UpdateIndices, Link12)
%linkClusters updates cluster IDs based on the linkages in Link12.
% This method will take the linkage information provided in Link12 and
% reassign ConnectID as appropriate to ensure linked localizations
% share the same ID.
%
% INPUTS:
%   ConnectID: Array of cluster IDs which indicate the connection between
%              localizations. (NLocalizations x 1 integer array)
%   MaxConnectID: Current maximum value of ConnectID.
%   UpdateIndices: SMD indices that are to be updated based on Link12.
%                  (NLocalizations x 1 integer array)
%   Link12: Array containing the trajectory link data as found by
%           solveLAP(). (NLocalizations x 1 integer array)
%
% OUTPUTS:
%   ConnectID: Updated connect IDs based on the input Link12.
%              (NLocalizations x 1 integer array)
%   MaxConnectID: Current maximum value of ConnectID.
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Initialize each localization as being a new cluster.
NLocalizations = numel(Link12) / 2;
ConnectID(UpdateIndices) = MaxConnectID + (1:NLocalizations).';

% Loop through the first NLocalizations entries of Link12 and assign new 
% connect IDs as needed, noting that the symmetry of the cost matrix makes
% the last NLocalizations entries redundant.
ConnectIDCurrent = ConnectID(UpdateIndices);
ConnectIDCopy = ConnectIDCurrent;
for ii = 1:NLocalizations
    % Skip to the next iteration if this link indicates a birth.
    if (Link12(ii) > NLocalizations)
        continue
    end
    
    % Determine the new ID which this cluster should be given.
    NewID = min([ConnectIDCurrent(ii), ConnectIDCurrent(Link12(ii)), ...
        ConnectIDCopy(ii), ConnectIDCopy(Link12(ii))]);
    
    % Update the connect IDs to reflect the connections.
    ConnectIDCurrent([ii, Link12(ii)]) = NewID;
end

% Store the updated connect IDs, ensuring we don't use a previously used
% ID.
ConnectID(UpdateIndices) = ConnectIDCurrent + MaxConnectID;
MaxConnectID = MaxConnectID + NLocalizations;


end