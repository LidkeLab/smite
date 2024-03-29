function [TR] = convertSMDToTR(SMD)
%convertSMDToTR converts an SMD into a TR structure.
% This method takes an SMD structure (see smi_core.SingleMoleculeData) and
% converts it into a TR structure.
% 
% INPUTS:
%   SMD: Single Molecule Data structure with a properly populated field
%        'ConnectID'.
% 
% OUTPUTS:
%   TR: Tracking Results.  The TR structure is a structure array containing
%       several fields present in the SMD structure but reorganized such
%       that each element of TR, e.g., TR(7), corresponds to a single
%       unique trajectory.
% 
% CITATION:

% Created by:
%   Hanieh Mazloom-Farsibaf (Lidke Lab, 2018)
%   Rewritten in smite, David J. Schodt (Lidke Lab, 2021)


% Create an empty TR structure, ending after this piece of code if no SMD
% was input (it might sometimes be useful to produce an empty TR
% structure).
TR = smi_core.TrackingResults.createTR();
if (~exist('SMD', 'var') || isempty(SMD))
    return
end

% Loop through each trajectory in SMD and place it in the output TR.
% NOTE: Some fields will be the same for all trajectories and are added
%       outside of the loop.
UniqueTrajIDs = unique(SMD.ConnectID);
NLocalizations = numel(SMD.FrameNum);
SMDFields = fieldnames(SMD);
NSMDFields = numel(SMDFields);
for ii = numel(UniqueTrajIDs):-1:1
    % Create an index array corresponding only to the current trajectory.
    CurrentTrajIndices = find(SMD.ConnectID == UniqueTrajIDs(ii));
    
    % Loop through each relevant SMD field and place it in the TR.
    for ff = 1:NSMDFields
        if (numel(SMD.(SMDFields{ff})) == NLocalizations)
            TR(ii, 1).(SMDFields{ff}) = ...
                SMD.(SMDFields{ff})(CurrentTrajIndices);
        else
            TR(ii, 1).(SMDFields{ff}) = SMD.(SMDFields{ff});
        end
    end
    
    % Store some other fields in the TR.
    TR(ii, 1).ConnectID = SMD.ConnectID(CurrentTrajIndices(1));
    TR(ii, 1).IndSMD = CurrentTrajIndices;
end


end