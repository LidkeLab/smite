function [TR] = threshTrajLength(TR, MinTrackLength)
%threshTrajLength removes trajectories smaller than a minimum track length.
% Given a TR and a MinTrackLength, this method will remove trajectories
% which are shorter than MinTrackLength from the TR structure, returning
% the updated TR structure with only sufficiently long trajectories.
%
% INPUTS:
%   TR: Tracking results structure containing information about single 
%       particle localizations and their trajectories.
%   MinTrackLength: Minimum trajectory length that will be kept in the
%                   output TR structure.
%   
% OUTPUTS:
%   TR: TR structure with short trajectories now removed. 
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2018)


% Remove trajectories with too few observations.
TrajLengths = smi_core.TrackingResults.computeTrajLengths(TR);
TR(TrajLengths < MinTrackLength) = [];


end