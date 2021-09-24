function [Trajectories, PeriodicityMap, ConnectionMap] = ...
    applyBoundaryCondition(...
    Trajectories, BoundaryCondition, FrameSize, ConnectionMap)
%applyBoundaryCondition applies a boundary condition to Trajectories.
% This method applies the specified 'BoundaryCondition' to the trajectories
% in 'Trajectories'.  For example, for a periodic boundary condition, the
% trajectories will be modified to ensure that localizations outside of
% the frame are periodically wrapped to the opposite side of the frame.
%
% NOTE: The method enforcePeriodicBoundary() can be used after this
%       method is applied to ensure trajectories that interact with a
%       periodic boundary are split into two distinct trajectories (i.e.,
%       two different rows in the array 'Trajectories').
%
% INPUTS:
%   Trajectories: Array of 2D trajectories. (NTrajx1x2)
%   BoundaryCondition: Boundary condition to be applied. (char array)
%                      'free': no changes are made at boundaries.
%                      'reflecting': trajectories bounce off of boundaries.
%                      'periodic': trajectories are periodically wrapped to
%                                  the opposite side of the frame.
%   FrameSize: Size of the frame defining the boundaries.
%              (pixels)(2x1 array)
%   ConnectionMap: Index mapping indicating oligomerization between
%                  trajectories. (see, e.g.,
%                  smi_sim.SimSPT.simOligomers() for usage)
%                  (NTrajx1)
% OUTPUTS:
%   Trajectories: Input 'Trajectories' updated by the boundary conditions.
%   PeriodicityMap: Logical array indicating trajectories for which a
%                   periodic boundary condition was applied. (NTrajx1)
%   ConnectionMap: Input 'ConnectionMap' after breaking up oligomers that
%                  hit a boundary.

% Created by:
%   David J. Schodt (Lidke Lab, 2021) with boundary condition usage
%       modified from the script InteractingParticles.m (unknown author(s)
%       in the Lidke Lab)


% Apply the boundary condition.
switch lower(BoundaryCondition)
    case 'periodic'
        TrajModFrameSize = ...
            [mod(Trajectories(:, 1, 1), FrameSize(1)), ...
            mod(Trajectories(:, 1, 2), FrameSize(2))];
        PeriodicityMap = any(TrajModFrameSize ...
            ~= reshape(Trajectories(:, 1, :), [], 2), 2);
        ConnectionMap = ConnectionMap .* ~PeriodicityMap;
        Trajectories(:, 1, 1) = TrajModFrameSize(:, 1);
        Trajectories(:, 1, 2) = TrajModFrameSize(:, 2);
    case 'reflecting'
        Trajectories(:, 1, 1) = Trajectories(:, 1, 1) ...
            + 2*((FrameSize(1)-Trajectories(:, 1, 1)) ...
            .*(Trajectories(:, 1, 1)>FrameSize(1)));
        Trajectories(:, 1, 2) = Trajectories(:, 1, 2) ...
            + 2*((FrameSize(2)-Trajectories(:, 1, 2)) ...
            .*(Trajectories(:, 1, 2)>FrameSize(2)));
        Trajectories(:, 1, :) = Trajectories(:, 1, :) ...
            - 2*Trajectories(:, 1, :).*(Trajectories(:, 1, :)<0);
end


end