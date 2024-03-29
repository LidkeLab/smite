function [TrajStruct] = simTrajectories(SimParams)
%simTrajectories simulates trajectories with oligomerization.
% This method simulates Brownian trajectories which may interact with one
% another to form dimers/higher order oligomers.
%
% INPUTS:
%   SimParams: Structure of simulation parameters.
%              (see smi_sim.SimSPT.defineDefaultParams())
%
% OUTPUTS:
%   TrajStruct: Structure containing information about the simulated
%               trajectories.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Ensure that 'SimParams' is complete, i.e., has all parameters.
if (~exist('SimParams', 'var') || isempty(SimParams))
    SimParams = smi_sim.SimSPT.defineDefaultParams();
else
    SimParams = smi_helpers.padStruct(SimParams, ...
        smi_sim.SimSPT.defineDefaultParams());
end

% Scatter the targets uniformly throughout the desired frame.
NTraj = poissrnd(SimParams.ParticleDensity * prod(SimParams.FrameSize));
InitialPositions = rand(NTraj, 2) .* SimParams.FrameSize;

% Apply the initial density mask.
InitialPositions = smi_sim.SimSPT.applyCoordMask(InitialPositions, ...
    SimParams.InitialDensityMask, SimParams.FrameSize);

% Simulate the trajectories and, if needed, oligomerization between them.
TrajStruct = smi_sim.SimSPT.simTrajBrownian(InitialPositions, SimParams);


end