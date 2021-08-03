function [DimerDurations] = ...
    computeDimerDurations(StateSequence, FrameNum, DimerState)
%computeDimerDurations computes the length of dimer events in StateSequence
% This method will look for changes in the values of StateSequence to
% determine when and for how long DimerState was present.  
% NOTE: If StateSequence(end) == DimerState, DimerDurations will treat that
%       as the end of a dimerization event.
%
% INPUTS:
%   StateSequence: An estimated sequence of states, e.g. [2, 1, 1, 2],
%                  which would correspond to an event length of 2.
%   FrameNum: Frame number corresponding to the observations in
%             StateSequence.
%   DimerState: Numeric value corresponding to the dimer state in
%               StateSequence.  If StateSequence = [1, 1, 1, 2, 2, 1, 2],
%               DimerState == 2 leads to DimerDurations = 2, and for
%               DimerState == 1 leads to DimerDurations = [3; 1].
%               (default = 1)
% 
% OUTPUTS:
%   DimerDurations: An array containing the observed durations of each
%                   dimerization event in StateSequence.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Set default parameters if needed.
if (~exist('DimerState', 'var') || isempty(DimerState))
    DimerState = 1;
end

% Ensure StateSequence is a column vector.
if (size(StateSequence, 1) < size(StateSequence, 2))
    StateSequence = StateSequence.';
end

% Create a boolean array to indicate which elements in StateSequence are
% equal to DimerState.
DimerEventBool = (StateSequence == DimerState);

% Find the starting/ending indices of each dimer event in StateSequence.
EventChanges = [DimerEventBool(1); diff(DimerEventBool)];
StartIndices = find(EventChanges == 1);
EndIndices = find(EventChanges == -1) - 1;

% If the start and end indices don't match in size, the last dimer event
% continued until the last observation.
if (numel(StartIndices) ~= numel(EndIndices))
    EndIndices = [EndIndices; numel(StateSequence)];
end

% Compute the duration of each event.
% NOTE: The +1 and the -1 above cancel each other out, but I like leaving
%       them because it helps clarify what the code is doing 
%       (at least for me)
DimerDurations = FrameNum(EndIndices) - FrameNum(StartIndices) + 1;


end