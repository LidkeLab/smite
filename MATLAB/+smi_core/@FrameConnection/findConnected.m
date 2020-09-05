function [SMDIndex] = findConnected(SMR, SMD, ID)
% findConnected finds localizations that were frame connected.
%
% Let N = numel(SMD.X) and M = numel(SMD_combined.X).  Then,
% i = find(SMD.ConnectID == id) where id = 1, ..., M produces the indices of
%    the SMD localizations that have been combined for each id.  Note id is the
%    number of a cluster in the internal frame connection algorithm and has
%    nothing to do with array indices.  Each localization in SMD_combined
%    corresponds to one cluster number, while one or more localizations in SMD
%    correspond to the same cluster number.
% j = find(SMD_combined.CombinedID == id) where id = 1, ..., M produces the
%    index of the SMD_combined localization that corresponds to the above
%    cluster id.
% Therefore, using the same id in both invocations produces the indices of the
% SMD localizations (i) that have been combined [via ConnectID] into a combined
% localization (j) [via CombinedID].
%
% INPUTS:
%    SMR:   SMR or SMD_combined structure from a frame-connected SMA_SR result
%    SMD:   SMD structure that contains threshold info, etc.
%    ID:    an index for a localization in a frame-connected SMR result
% OUTPUTS:
%    SMDIndex:  An array of indices for the SMD structure for the localizations
%               that were connected to form the localization ID in SMR. 

% Created by
%    Michael Wester (LidkeLab 2020)

   clusterID = SMR.ConnectID(ID);
   SMDIndex = find(SMD.ConnectID == clusterID);

end