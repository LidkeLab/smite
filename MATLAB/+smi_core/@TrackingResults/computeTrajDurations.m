function [Durations] = computeTrajDurations(TR)
%computeTrajDurations computes the duration of trajectories in TR.
% This method computes the duration of all trajectories in the provided TR
% structure.
% 
% INPUTS:
%   TR: Tracking Results structure.
% 
% OUTPUTS:
%   Durations: Duration of each trajectory in TR.  The indexing is matched
%              to TR, e.g., trajectory TR(m) had a duration of Durations(m)
%              frames. (frames)
% 
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Compute the trajectory durations.
Durations = cellfun(@(X) max(X) - min(X) + 1, {TR.FrameNum});


end