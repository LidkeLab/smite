function [Fidelity] = computeTrajFidelity(TR)
%computeTrajFidelity computes the trajectory length divided by duration.
% This method computes the length (number of observations) in each
% trajectory in TR divided by the duration (max-min frame + 1).
% 
% INPUTS:
%   TR: Tracking Results structure.
% 
% OUTPUTS:
%   Fidelity: Fidelity (length/duration) for each trajectory in TR.
% 
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Loop through trajectories and compute the trajectory fidelities.
Fidelity = NaN(size(TR));
for nn = 1:numel(TR)
    Fidelity(nn) = smi_core.TrackingResults.computeTrajLengths(TR(nn)) ...
        / smi_core.TrackingResults.computeTrajDurations(TR(nn));
end


end