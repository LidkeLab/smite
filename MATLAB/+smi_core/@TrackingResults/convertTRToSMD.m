function [SMD] = convertTRToSMD(TR)
%convertTRToSMD converts a TR back into an SMD structure.
% This method acts as the inverse of convertSMDToTR().
% 
% INPUTS:
%   TR: Tracking Results.  The TR structure is a structure array containing
%       several fields present in the SMD structure but reorganized such
%       that each element of TR, e.g., TR(7), corresponds to a single
%       unique trajectory.
% 
% OUTPUTS:
%   SMD: Single Molecule Data structure with a properly populated field
%        'ConnectID'.
% 
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Create an empty SMD structure, ending after this piece of code if no TR
% was input (it might sometimes be useful to produce an empty SMD
% structure).
SMD = smi_core.SingleMoleculeData.createSMD();
if (~exist('TR', 'var') || isempty(TR))
    return
end

% Loop through each trajectory in TR and place its localizations in the
% output SMD.
for ii = 1:length(TR)
    SMD = smi_core.SingleMoleculeData.catSMD(SMD, TR(ii), false, false);
end

% Fix some fields that don't exactly match our expectations.
SMD.IndSMD = {};
SMD.ConnectID = ones(size(SMD.FrameNum));
NLoc = smi_core.TrackingResults.computeTrajLengths(TR);
NLocCumulative = [0; cumsum(NLoc)];
for ii = 1:length(TR)
    IndArray = (1:NLoc(ii)) + NLocCumulative(ii);
    SMD.ConnectID(IndArray) = TR(ii).ConnectID * ones(NLoc(ii), 1);
end


end