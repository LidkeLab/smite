function [SMD] = connectTrajGC(SMD, Link12)
%connectTrajGC connects gaps in trajectories in an SMD structure.
% Given an Mx1 vector Link12 containing gap closing data found by
% smi.SPT.solveLAP(), this method will connect (as appropriate)
% trajectories in a global sense.
%
% NOTE: This method assumes that min(SMD.ConnectID)==1 .
%
% INPUTS:
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)
%        containing the localizations that we wish to stitch into
%        trajectories.
%   Link12: Vector containing the gap closing data as found by
%           smi.SPT.solveLAP(), where M is twice the number of trajectories
%           which have been considered by the gap closing LAP. 
%           Link12 is defined such that, e.g., if Link12(7) = 3, 
%           the end of the third trajectory should be linked to the start 
%           of the seventh trajectory. (Mx1)
%
% OUTPUTS:
%   SMD: Input SMD structure which contains the modified ConnectID field 
%        corresponding to newly linked localizations.
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke lab, 2018)
%   Revised by Hanieh Mazloom-Farsibaf (Lidke lab, 2019)
%   Reorganized with minor revisions, David J. Schodt (Lidke Lab, 2020)


% Keep original ConnectID before Gap Closing (=SMD.PreGCConnectID) and
% define the SMD.ConnectID for final results
SMD.PreGCConnectID = SMD.ConnectID;

% Determine the maximum value of ConnectID for later use.
MaxTrajID = max(SMD.ConnectID);

% Loop through all of the potential gaps, assuming that the ConnectID
% consists of consecutive integers up until this point.
for ee = 1:MaxTrajID
    % If a gap has been closed, trajectories with ConnectID's Link12(ee)
    % and ee (which were previously thought to be distinct trajectories)
    % should be linked, and thus their ConnectID's should be set to the
    % same value. The linked trajectories will be assigned a ConnectID
    % that is the minimum of the ConnectID of the trajectory end being
    % considered (which is just ee), the ConnectID of the trajectory
    % beginning it was linked to (Link12(ee)), and the ID which was 
    % previously re-assigned to another segment associated with the ee-th
    % (or Link12(ee)-th) trajectory.
    
    % If the end of trajectory ee was linked to the start of
    % another trajectory (and not to a death), we will assign each of them
    % the same trajectory ID.
    if (Link12(ee) <= MaxTrajID)
        % Create a boolean array indicating which trajectory ID's we will
        % be updating.
        UpdateBoolean = (SMD.ConnectID == max(ee, Link12(ee)));
                
        % If the ee-th trajectory ID was already re-assigned, we need to
        % account for that as well.
        ReassignedCurrentID = min(SMD.ConnectID(SMD.PreGCConnectID == ee));
        
        % If linking to a trajectory which itself had already been given a
        % new trajectory ID (i.e. the trajectory with id Link12(ee) has
        % been assigned a new ConnectID ~= Link12(ee)), we also want to
        % consider assigning the current trajectory the updated ConnectID
        % associated with Link12(ee).
        ReassigedLinkID = ...
            min(SMD.ConnectID(SMD.PreGCConnectID == Link12(ee)));
        
        % Re-assign the appropriate set of trajectories to their new
        % ConnectID.
        SMD.ConnectID(UpdateBoolean) = ...
            ones(sum(UpdateBoolean), 1, 'int32') ...
            * min([ee, Link12(ee), ReassignedCurrentID, ReassigedLinkID]);
    end
end

% Ensure SMD.ConnectID consists of the set of integers 1:NTraj without
% skipping any integers.
SMD.ConnectID = smi_helpers.compressToRange(SMD.ConnectID);


end