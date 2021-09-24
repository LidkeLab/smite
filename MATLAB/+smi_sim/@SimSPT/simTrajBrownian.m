function [TrajStruct] = simTrajBrownian(InitialPositions, SimParams)
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
%   TrajStruct: Structure containing information about the simulated
%               trajectories.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Make local copies of some parameters (can improve speed for large
% simulations).
FrameSize = SimParams.FrameSize;
D = SimParams.D;
BoundaryCondition = SimParams.BoundaryCondition;
InteractionDistance = SimParams.InteractionDistance;
InteractionProb = SimParams.InteractionProb;
RestrictToDimers = SimParams.RestrictToDimers;

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
KDisconnectSub = SimParams.KDisconnect / SimParams.SubframeDensity;
NSubframes = SimParams.NFrames * SimParams.SubframeDensity;
Trajectories = zeros(NTraj, NSubframes, 2, 'double');
Trajectories(:, 1, :) = InitialPositions;
PeriodicityMapT = zeros(NTraj, NSubframes);
ConnectionMap = zeros(NTraj, 1);
ConnectionMapT = zeros(NTraj, NSubframes);
IsOligoSim = (SimParams.InteractionProb ...
    && ~isinf(SimParams.InteractionDistance));
for ff = 2:NSubframes
    % Sample the proposed trajectory updates from the Normal
    % distribution (Brownian motion).
    TrajectoryUpdates = sqrt(2*DSub) .* randn(NTraj, 1, 2);
    
    % Simulate oligomerization.
    if IsOligoSim
        [Trajectories(:, ff, :), ConnectionMap] = ...
            smi_sim.SimSPT.simOligomers(...
            Trajectories(:, ff-1, :), TrajectoryUpdates, NTraj, ...
            ConnectionMap, InteractionDistance, InteractionProb, ...
            KDisconnectSub, RestrictToDimers);
    end
    
    % Apply the boundary conditions.
    [Trajectories(:, ff, :), ...
        PeriodicityMapT(:, ff), ConnectionMap] = ...
        smi_sim.SimSPT.applyBoundaryCondition(Trajectories(:, ff, :), ...
        BoundaryCondition, FrameSize, ConnectionMap);
    ConnectionMapT(:, ff) = ConnectionMap;
end

% Break up trajectories that experienced a periodic boundary (that is, each
% interaction with the periodic boundary creates a new trajectory starting
% from that interaction).
[Trajectories, ConnectionMapT, IsOn, TrajMap] = ...
    smi_sim.SimSPT.enforcePeriodicBoundary(...
    Trajectories, PeriodicityMapT, ConnectionMapT);
TrajStruct.IsOn = IsOn;
TrajStruct.D = DSub(TrajMap);
DataSize = size(IsOn);
TrajStruct.Photons = zeros(DataSize);
TrajStruct.Photons_SE = inf(DataSize);
TrajStruct.Bg = zeros(DataSize);
TrajStruct.Bg_SE = inf(DataSize);
TrajStruct.ConnectionMapT = ConnectionMapT;
TrajStruct.Trajectories = Trajectories;
TrajStruct.Trajectories_SE = zeros([DataSize, 2]);


end