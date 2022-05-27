function [NObservations] = computeTrajLengths(TR)
%computeTrajLengths computes the length of trajectories in TR.
% This method computes the length (number of observations) in each
% trajectory in TR.
% 
% INPUTS:
%   TR: Tracking Results structure.
% 
% OUTPUTS:
%   NObservations: Length of each trajectory in TR.  The indexing is 
%                  matched to TR, e.g., trajectory TR(m) had a duration of
%                  NObservations(m) frames. (frames)
% 
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Compute the trajectory lengths.
NObservations = cellfun(@numel, {TR.FrameNum}.');


end