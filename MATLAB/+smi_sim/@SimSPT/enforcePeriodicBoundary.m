function [Trajectories, ConnectionMapT, IsOn, TrajMap] = ...
    enforcePeriodicBoundary(Trajectories, PeriodicityMapT, ConnectionMapT)
%enforcePeriodicBoundary breaks trajectories that hit a periodic boundary.
% If a trajectory in 'Trajectories' experienced a periodic boundary, this
% method will break that trajectory into two trajectories: one before and
% one after the periodic boundary (i.e., when the particle is translated to
% the other side of the frame, it is treated as a new trajectory).  This
% method will also break apart oligomer events when one of the constituents
% experiences a periodic boundary.
%
% INPUTS:
%   Trajectories: Array of 2D trajectories. (NTrajxNFramesx2)
%   PeriodicityMapT: Binary map indicating when 'Trajectories' experienced
%                    a periodic boundary (1 indicates a periodic boundary
%                    was encountered for that trajectory in that frame).
%                    (see, e.g., smi_sim.SimSPT.simOligoTrajBrownian() for
%                    usage)(NTrajxNFrames)
%   ConnectionMapT: Index mapping indicating oligomerization between
%                   trajectories. (see, e.g.,
%                   smi_sim.SimSPT.simOligoTrajBrownian() for usage)
%                   (NTrajxNFrames)
% OUTPUTS:
%   Trajectories: Input 'Trajectories' after breaking up the trajectories
%                that encountered the periodic boundaries.
%   ConnectionMapT: Input 'ConnectionMapT' after breaking up oligomers that
%                   hit a boundary.
%   IsOn: Flag indicating when the trajectories are "on" (in this case, the
%         post-birth and pre-death frames of a trajectory).  This is
%         returned to improve speed for large simulations, since checking
%         isnan(Trajectories) might be slower than just keeping tracking of
%         those entries as they're made. (NTrajxNFrames)
%   TrajMap: Set of indices mapping the rows of output 'Trajectories' to 
%            their original row in the input 'Trajectories' (e.g., if the
%            input TrajectoriesInput(5, :, :) was split by a periodic
%            boundary into TrajectoriesOutput(5, :, :) and 
%            TrajectoriesOutput(12, :, :), TrajMap(5) = 5 and
%            TrajMap(12) = 5.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('ConnectionMapT', 'var') || isempty(ConnectionMapT))
    ConnectionMapT = zeros(size(PeriodicityMapT));
end

% Enforce the periodic boundaries, setting elements of 'Trajectories' to
% NaN before new trajectories are birthed and after old ones die.
NTrajInitial = size(Trajectories, 1);
NFrames = size(Trajectories, 2);
CurrentTrajID = NTrajInitial;
IsOn = ones(size(PeriodicityMapT), 'logical');
TrajMap = (1:NTrajInitial).';
for dd = 1:NTrajInitial
    % Loop through any periodic boundary encounters and turn trajectories 
    % off in the appropriate frames.
    if any(PeriodicityMapT(dd, :))
        EventIndices = find(PeriodicityMapT(dd, :));
        NNewTraj = numel(EventIndices);
        NewTrajectories = repmat(Trajectories(dd, :, :), [NNewTraj, 1]);
        NewConnections = repmat(ConnectionMapT(dd, :), [NNewTraj, 1]);
        IsOn(dd, (EventIndices(1):NFrames)) = false;
        Trajectories(dd, (EventIndices(1):NFrames), :) = NaN;
        ConnectionMapT(dd, (EventIndices(1):NFrames)) = NaN;
        for nn = 1:(NNewTraj-1)
            % Update the trajectory ID counter.
            CurrentTrajID = CurrentTrajID + 1;
            
            % Ensure this particle is off before its birth and after
            % its death.
            IsOn(CurrentTrajID, 1:(EventIndices(nn)-1)) = false;
            IsOn(CurrentTrajID, EventIndices(nn+1):NFrames) = false;
            
            % Revise the trajectory positions and connections to
            % reflect the birth.
            NewTrajectories(nn, 1:(EventIndices(nn)-1), :) = NaN;
            NewTrajectories(nn, EventIndices(nn+1):NFrames, :) = NaN;
            NewConnections(nn, 1:(EventIndices(nn)-1)) = NaN;
            NewConnections(nn, EventIndices(nn+1):NFrames) = NaN;
        end
        CurrentTrajID = CurrentTrajID + 1;
        IsOn(CurrentTrajID, 1:(EventIndices(NNewTraj)-1)) = false;
        NewTrajectories(NNewTraj, 1:(EventIndices(NNewTraj)-1), :) = NaN;
        Trajectories = cat(1, Trajectories, NewTrajectories);
        TrajMap = [TrajMap; repelem(dd, NNewTraj, 1)];
        NewConnections(NNewTraj, 1:(EventIndices(NNewTraj)-1)) = NaN;
        ConnectionMapT = cat(1, ConnectionMapT, NewConnections);
    end
end


end