function [SMD] = connectTrajFF(SMD, Link12, FrameNumber)
%connectTrajFF connects frame-to-frame localizations into trajectories.
% Given an (M+N) x 1 vector containing localization link data found by
% smi.SPT.solveLAP(), this method will connect (as appropriate) the N
% localizations in frame FrameNumber to the M localizations in frame 
% FrameNumber+1.
%
% INPUTS:
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)
%        containing the localizations that we wish to stitch into
%        trajectories.
%   Link12: Array containing the trajectory link data as found 
%           by smi.SPT.solveLAP(), where M+N is the sum of the number of 
%           localizations in frame FrameNumber+1 (M) and the number of 
%           localizations in frame FrameNumber (N). 
%           For example, if Link12(4) = 1, we need to link the fourth
%           localization in frame FrameNumber+1 to the first localization 
%           in frame FrameNumber. (M+N x 1) 
%   FrameNumber: The frame whose localizations will be connected to the 
%                subsequent frame FrameNumber+1. (1x1) 
%
% OUTPUTS:
%   SMD: Input SMD structure which contains the modified ConnectID field 
%        corresponding to newly linked localizations.
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2018)


% Determine the total number of localizations in each of the two frames.
NLinks = numel(Link12);
NLocFirstFrame = sum(SMD.FrameNum == FrameNumber);
NLocSecondFrame = NLinks - NLocFirstFrame;

% Grab the maximum trajectory ID from the TD structure for later use.  If
% the field ConnectID is empty, initialize it first to an array of zeros
if isempty(SMD.ConnectID)
    SMD.ConnectID = zeros(numel(SMD.X), 1);
end
MaxTrajID = max(SMD.ConnectID);

% Find the set of indices corresponding to localizations in frames
% FrameNumber and FrameNumber+1 (frame 1 and frame 2).
FrameOneIndices = find(SMD.FrameNum == FrameNumber); 
FrameTwoIndices = find(SMD.FrameNum == (FrameNumber+1));

% Loop through each localization in frame FrameNumber and check if it has
% been associated with a ConnectID.  If it has not, set it's
% ConnectID to be the smallest possible non-existing ConnectID.
for nn = 1:NLocFirstFrame
    if ~SMD.ConnectID(FrameOneIndices(nn))
        MaxTrajID = MaxTrajID + 1; 
        SMD.ConnectID(FrameOneIndices(nn)) = MaxTrajID; 
    end
end

% Loop through the localizations in frame FrameNumber+1 and assign them a
% trajectory ID based on Link12.
% NOTE: We only care about the first NTrajSecondFrame elements of Link12
%       because the rest correspond to the "death" and "auxillary" blocks
%       of the cost matrix.
for nn = 1:NLocSecondFrame
    % Determine if Link12(nn) is linking this localization to a
    % localization in the first frame.  If it's not, assign to this
    % localization the smallest unique trajectory ID available.
    if (Link12(nn) > NLocFirstFrame)
        MaxTrajID = MaxTrajID + 1;
        SMD.ConnectID(FrameTwoIndices(nn)) = MaxTrajID;
    else
        SMD.ConnectID(FrameTwoIndices(nn)) = ...
            SMD.ConnectID(FrameOneIndices(Link12(nn)));
    end
end


end