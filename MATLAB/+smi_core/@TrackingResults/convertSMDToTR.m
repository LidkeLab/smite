function [TR] = convertSMDToTR(SMD, FileInfoStruct)
%convertSMDToTR converts an SMD into a TR structure.
% This method takes an SMD structure (see smi_core.SingleMoleculeData) and
% converts it into a TR structure.
% 
% INPUTS:
%   SMD: Single Molecule Data structure with a properly populated field
%        'ConnectID'.
%   FileInfoStruct: A structure array with fields FileDir and FileName
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


% Set defaults as needed.
if (~exist('FileInfoStruct', 'var') || isempty(FileInfoStruct))
    FileInfoStruct.FileDir = '';
    FileInfoStruct.FileName = '';
end

% Create an empty TR structure, ending after this piece of code if no SMD
% was input (it might sometimes be useful to produce an empty TR
% structure).
TR = smi_core.TrackingResults.createTR();
[TR.FileDir] = deal(FileInfoStruct.FileDir);
[TR.FileName] = deal(FileInfoStruct.FileName);
SMDFields = fieldnames(SMD);
if (~exist('SMD', 'var') || isempty(SMD) || isempty(SMDFields))
    return
end

% Loop through each trajectory in SMD and place it in the output TR.
% NOTE: Some fields will be the same for all trajectories and are added
%       outside of the loop.
UniqueTrajIDs = unique(SMD.ConnectID);
NLocalizations = numel(SMD.FrameNum);
TRFields = fieldnames(TR);
TRFields = TRFields(ismember(TRFields, SMDFields));
NTRFields = numel(TRFields);
for ii = numel(UniqueTrajIDs):-1:1
    % Create an index array corresponding only to the current trajectory.
    CurrentTrajIndices = find(SMD.ConnectID == UniqueTrajIDs(ii));
    
    % Loop through each relevant SMD field and place it in the TR.
    for ff = 1:NTRFields
        if (numel(SMD.(TRFields{ff})) == NLocalizations)
            TR(ii, 1).(TRFields{ff}) = SMD.(TRFields{ff})(CurrentTrajIndices);
        else
            TR(ii, 1).(TRFields{ff}) = SMD.(TRFields{ff});
        end
    end
    
    % Store some other fields in the TR.
    TR(ii, 1).ConnectID = SMD.ConnectID(CurrentTrajIndices(1));
    TR(ii, 1).IndSMD = CurrentTrajIndices;
end
[TR.FileDir] = deal(FileInfoStruct.FileDir);
[TR.FileName] = deal(FileInfoStruct.FileName);


end