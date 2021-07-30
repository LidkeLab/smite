function [TRArrayTrunc] = isolateCandidateTRArray(TRArray)
%isolateCandidateTRArray isolates dimer candidate pairs in TRArray.
% This method will loop through each trajectory in TRArray and isolate the
% portions of the trajectory considered candidates for dimer events.
% Specifically, the input TRArray should have a field DimerCandidateBool,
% and this method will simply loop through all fields and apply the
% DimerCandidateBool mask to the data.
%
% INPUTS:
%   TRArray: A structure array of TR structures, where the constituent TR
%            structures correspond to dimer candidate trajectories.  The 
%            first index corresponds to the "channel" of the trajectory and 
%            the second index corresponds to the pair number.  For example,
%            TRArray(1, j) will contain information about a trajectory from 
%            TR1 that was observed within MaxDimerDistance of
%            TRArray(2, j), a trajectory in TR2.
%
% OUTPUTS:
%   TRArrayTrunc: The input structure TRArray with all fields truncated to
%                 only show the portions of the trajectories that were
%                 marked as dimer candidates.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Generate a list of fieldnames in the TRArray.
TRArrayFields = fieldnames(TRArray);

% Loop through all fields and apply the DimerCandidateBool mask where
% relevant.
TRArrayTrunc = TRArray;
for ii = 1:size(TRArray, 2)
    DimerCandidateBool1 = logical(TRArray(1, ii).DimerCandidateBool);
    DimerCandidateBool2 = logical(TRArray(2, ii).DimerCandidateBool);
    for jj = 1:numel(TRArrayFields)
        % Determine how many datapoints exist in each of the two
        % trajectories in this pairing.
        NDataPoints1 = numel(TRArray(1, ii).FrameNum);
        NDataPoints2 = numel(TRArray(2, ii).FrameNum);
        
        % Apply the DimerCandidateBool mask to relevant fields.
        if (numel(TRArray(1, ii).(TRArrayFields{jj})) == NDataPoints1)
            TRArrayTrunc(1, ii).(TRArrayFields{jj}) = ...
                TRArray(1, ii).(TRArrayFields{jj})(DimerCandidateBool1);
        end
        if (numel(TRArray(2, ii).(TRArrayFields{jj})) == NDataPoints2)
            TRArrayTrunc(2, ii).(TRArrayFields{jj}) = ...
                TRArray(2, ii).(TRArrayFields{jj})(DimerCandidateBool2);
        end
    end
end


end