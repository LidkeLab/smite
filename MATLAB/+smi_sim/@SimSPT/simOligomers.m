function [Trajectories, ConnectionMap] = simOligomers(...
    Trajectories, TrajectoryUpdates, NTraj, ...
    ConnectionMap, InteractionDistance, InteractionProb, KDisconnect, ...
    RestrictToDimers)
%simOligomers simulates oligomerization between trajectories.
% This method updates 'Trajectories' based on their oligomerization state
% and also creates and destroys oligomers as appropriate.
%
% INPUTS:
%   Trajectories: Array of 2D trajectories. (NTrajx1x2)
%   TrajectoryUpdates: Proposed steps to be made by the 'Trajectories'.
%                      (NTrajx1x2)
%   NTraj: Number of trajectories, passed to improve speed instead of
%          recomputing each time this method is called. (scalar)
%   ConnectionMap: Index mapping indicating oligomerization between
%                  trajectories. For example, if ConnectionMap(mm)==nn,
%                  Trajectories(mm, :, :) is in an oligomer with
%                  Trajectories(nn, :, :). (NTrajx1)
%   InteractionDistance: Distance at or below which trajectories have a
%                        chance to interact with one another. (scalar)
%   InteractionProb: Probability that two trajectories interact when at or
%                    below 'InteractionDistance'. (scalar)
%   KDisconnect: Rate parameter for oligomer disconnection. (scalar)
%   RestrictToDimers: Flag indicating that only dimers are allowed (i.e.,
%                     no higher order oligomers are allowed between
%                     trajectories). (Boolean scalar)
%
% OUTPUTS:
%   Trajectories: Input 'Trajectories' updated based on oligomerization
%                 state and 'TrajectoryUpdates'.
%   ConnectionMap: Input 'ConnectionMap' after updating to reflect new or
%                  broken oligomers.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Force oligomers to move together.
for mm = 1:NTraj
    if ConnectionMap(mm)
        TrajectoryUpdates(mm, :) = TrajectoryUpdates(ConnectionMap(mm), :);
    end
end

% Update the particle positions.
Trajectories = Trajectories + TrajectoryUpdates;

% Update the connection map.
OldConnectionMap = ConnectionMap;
for mm = 1:NTraj
    for nn = (mm+1):NTraj
        % Compute the distance between the two particles currently
        % being compared.
        DeltaX = Trajectories(mm, 1, 1) - Trajectories(nn, 1, 1);
        DeltaY = Trajectories(mm, 1, 2) - Trajectories(nn, 1, 2);
        ParticleSeparation = sqrt(DeltaX^2 + DeltaY^2);
        
        % Determine if these particles should be connected.
        if ((ParticleSeparation<=InteractionDistance) ...
                && (rand()<InteractionProb))
            % If these particles weren't connected in the previous
            % frame, then we will connect them now.
            if ~(RestrictToDimers || ConnectionMap(nn))
                % Oligomers of higher order than 2 are allowed.
                ConnectionMap(nn) = mm;
            elseif ~(ConnectionMap(nn) || ConnectionMap(mm))
                ConnectionMap(nn) = mm;
                ConnectionMap(mm) = nn;
            else
                continue
            end
            
            % Force these two particles to be exactly
            % InteractionDistance apart from one another.
            DistanceDiff = InteractionDistance - ParticleSeparation;
            Theta = atan(DeltaY / DeltaX);
            XShift = sign(DeltaX) * abs(DistanceDiff*cos(Theta)) / 2;
            YShift = sign(DeltaY) * abs(DistanceDiff*sin(Theta)) / 2;
            Trajectories(nn, 1, 1) = Trajectories(nn, 1, 1) - XShift;
            Trajectories(nn, 1, 2) = Trajectories(nn, 1, 2) - YShift;
            Trajectories(mm, 1, 1) = Trajectories(mm, 1, 1) + XShift;
            Trajectories(mm, 1, 2) = Trajectories(mm, 1, 2) + YShift;
        end
    end
end

% Randomly disconnect particle connections.
for nn = 1:NTraj
    if OldConnectionMap(nn)
        % Set the previous connection map for this connection to 0
        % (that way we can skip it when reached by the outer loop).
        OldConnectionMap(OldConnectionMap(nn)) = 0;
        
        % Randomly disconnect the particle connection by comparing
        % rand() to the probability that these particles will
        % disconnect.
        if (rand() < (1-exp(-KDisconnect)))
            % Ensure that the connection is broken for the next
            % frame.
            ConnectionMap(ConnectionMap(nn)) = 0;
            ConnectionMap(nn) = 0;
        end
    end
end


end