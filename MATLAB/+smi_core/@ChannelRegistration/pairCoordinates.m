function [PairMap12, PairMap21] = pairCoordinates(Coords1, Coords2, ...
    SeparationThreshold)
%pairCoordinates pairs sets of coordinates with each other.
% This method will pair the coordinates in Coords1 with the coordinates in
% Coords2 based on their distance (i.e., this finds the nearest neighbor of
% each of Coords1 in Coords2), ensuring the Threshold is note exceeded and
% that no duplicate pairings are made (each of Coords1 is only paired to
% one of Coords2, and vice versa).
%
% INPUTS:
%   Coords1: A set of x,y coordinates (numeric array, Mx2)
%   Coords2: A set of x,y coordinates to be paired with Coords1
%            (numeric array, Nx2)
%   SeparationThreshold: Maximum distance between paired coordinates.
%                        (Pixels)(positive scalar)(Default = inf)
% 
% OUTPUTS:
%   PairMap12: Array specifying how Coords1 and Coords2 were paired.  For
%              example, if PairMap12(7) = 5, Coords1(7, :) and 
%              Coords2(5, :) were paired together.  A NaN marker indicates 
%              a pairing was not made for that index.
%              (Px1 array, P<=min(M, N))
%   PairMap21: Array specifying how Coords1 and Coords2 were paired.  For
%              example, if PairMap21(3) = 12, Coords2(3, :) and 
%              Coords1(12, :) were paired together.  A NaN marker indicates 
%              a pairing was not made for that index.
%              (Px1 array, P<=min(M, N))

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Set defaults if needed.
if (~exist('SeparationThreshold', 'var') || isempty(SeparationThreshold))
    SeparationThreshold = inf;
end

% Loop through Coords1 and search for the nearest neighbor.
NCoords1 = size(Coords1, 1);
PairMap12 = NaN(NCoords1, 1);
PairMap21 = PairMap12;
for ii = 1:NCoords1
    % Find the nearest coordinate in Coords2 to Coords1(ii, :).
    Separation = sqrt((Coords1(ii, 1)-Coords2(:, 1)).^2 ...
        + (Coords1(ii, 2)-Coords2(:, 2)).^2);
    [MinSeparation, NearestNeighborIndex] = min(Separation);
    if ((MinSeparation<=SeparationThreshold) ...
            && ~ismember(NearestNeighborIndex, PairMap12))
        PairMap12(ii) = NearestNeighborIndex;
        PairMap21(NearestNeighborIndex) = ii;
    end
end


end