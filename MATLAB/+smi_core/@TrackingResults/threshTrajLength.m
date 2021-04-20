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


% Loop through each trajectory in the TR structure. 
% NOTE: It's shorter to do this with cellfun() (find TrackLength for all 
%       trajectories at once then threshold), but that's several times
%       slower than doing this loop below.
IndicesToKeep = []; % initialize our array of traj. indices we wish to keep
NTraj = numel(TR); % number of trajectories in the TR structure
for ii = 1:NTraj 
    % Determine the length of the ii-th trajectory in TR (number of frames
    % it exists in).
    TrackLength = numel(TR(ii).FrameNum); 
    
    % If the trajectory is sufficiently long, save it's index for later.
    % Otherwise, do nothing.
    if (TrackLength >= MinTrackLength)
    	IndicesToKeep = [IndicesToKeep; ii];
    end
end

% Isolate the trajectories that we wish to keep from the TR structure.
TR = TR(IndicesToKeep);


end