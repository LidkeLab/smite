function [ConnectID] = validifyConnectID(ConnectID)
%validifyConnectID ensures ConnectID consists of the numbers 1:NConnectID
% This method will convert an input ConnectID to the set of numbers
% 1:numel(ConnectID).  For example, if 
% ConnectID = [1, 2, 2, 4, 5, 5, 5, 7], the output will be 
% [1, 2, 2, 3, 4, 4, 4, 5].
%
% INPUT:
%   ConnectID: An array of integers.
%
% OUTPUT:
%   ConnectID: Array of integers with the same equality pattern as the
%              input, but consisting only of the integers
%              [1, numel(ConnectID)]

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Ensure ConnectID consists of the set of integers 1:numel(ConnectID)
% without skipping any integers.
ConnectIDList = unique(ConnectID);
NIDs = numel(ConnectIDList);
for nn = 1:NIDs
    % Set all entries with ConnectID==ConnectIDList(nn) equal to nn.
    ConnectID(ConnectID==ConnectIDList(nn)) = nn;
end


end