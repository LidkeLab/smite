function [TrajectoryStruct] = simTrajBrownian(InitialPositions, SimParams)
%simTrajBrownian simulates Brownian trajectories.
% This method simulates Brownian random motion of the targets at the
% specified initial conditions.
%
% INPUTS:
%   InitialPositions: Initial positions of the random walkers 
%                     (NTrajx2 array)
%   SimParams: Structure of simulation parameters (see
%              smi_sim.SimSPT.defineDefaultParams())
%
% OUTPUTS:
%   TrajectoryStruct: Structure containing information about the simulated
%                     trajectories.

% Created by:
%   David J. Schodt (Lidke Lab, 2021) with boundary condition usage
%       modified from the script InteractingParticles.m (unknown author(s)
%       in the Lidke Lab)


% Make local copies of some parameters (can improve speed for large
% simulations).
FrameSize = SimParams.FrameSize;
D = SimParams.D;
BoundaryCondition = SimParams.BoundaryCondition;

% Prepare the diffusion coefficients.
NTraj = size(InitialPositions, 1);
NDiffusionCoefficients = numel(D);
if (NDiffusionCoefficients > NTraj)
    D = D(1:NTraj);
else
    D = D(randi(NDiffusionCoefficients, [NTraj, 1]));
end
if isrow(D)
    D = D.';
end

% Simulate the Brownian trajectories.
DSub = D / SimParams.SubframeDensity;
NSubframes = SimParams.NFrames * SimParams.SubframeDensity;
Trajectories = zeros(NTraj, NSubframes, 2, 'double');
Trajectories(:, 1, :) = InitialPositions;
PeriodicityMapT = zeros(NTraj, NSubframes);
ConnectionMapT = zeros(NTraj, NSubframes);
for ff = 2:NSubframes
    % Update each trajectory by sampling the X, Y displacements from a
    % normal distribution.
    TrajectoryUpdates = sqrt(2*DSub) .* randn(NTraj, 1, 2);
    
    % Update the particle positions.
    Trajectories(:, ff, :) = Trajectories(:, ff-1, :) ...
        + TrajectoryUpdates;
    
    % Apply the boundary conditions.
    switch BoundaryCondition
        case 'Periodic'
            TrajModFrameSize = ...
                [mod(Trajectories(:, ff, 1), FrameSize(1)), ...
                mod(Trajectories(:, ff, 2), FrameSize(2))];
            PeriodicityMapT(:, ff) = any(TrajModFrameSize ...
                ~= squeeze(Trajectories(:, ff, :)), 2);
            Trajectories(:, ff, 1) = TrajModFrameSize(:, 1);
            Trajectories(:, ff, 2) = TrajModFrameSize(:, 2);
        case 'Reflecting'
            Trajectories(:, ff, 1) = Trajectories(:, ff, 1) ...
                + 2*((FrameSize(1)-Trajectories(:, ff, 1)) ...
                .*(Trajectories(:, ff, 1)>FrameSize(1)));
            Trajectories(:, ff, 2) = Trajectories(:, ff, 2) ...
                + 2*((FrameSize(2)-Trajectories(:, ff, 2)) ...
                .*(Trajectories(:, ff, 2)>FrameSize(2)));
            Trajectories(:, ff, :) = Trajectories(:, ff, :) ...
                - 2*Trajectories(:, ff, :).*(Trajectories(:, ff, :)<0);
    end
end

% Break up trajectories that experienced a periodic boundary (that is, each
% interaction with the periodic boundary creates a new trajectory starting
% from that interaction).
[Trajectories, ConnectionMapT, IsOn, TrajMap] = ...
    smi_sim.SimSPT.enforcePeriodicBoundary(...
    Trajectories, PeriodicityMapT, ConnectionMapT);
TrajectoryStruct.IsOn = IsOn;
TrajectoryStruct.DSub = DSub(TrajMap);
TrajectoryStruct.PhotonsSub = [];
TrajectoryStruct.ConnectionMapT = ConnectionMapT;
TrajectoryStruct.Trajectories = Trajectories;


end