function [SubPixelOffsets, SMDROIs, SMDStats] = ...
    estimateLocalCoordShifts(SMD1, SMD2, SubROISize)
%estimateLocalCoordShifts estimates local shifts between two SMDs.
%
% INPUT:
%   SMD1: Single Molecule Data structure containing localizations that will
%         be compared to SMD2 localizations.
%   SMD2: Single Molecule Data structure containing localizations that will
%         be compared to SMD1 localizations.
%   SubROISize: The size of local regions in which the shift will be
%               computed, ideally evenly divides [m, n].
%               (Pixels)(2x1 array)(Default = [SMD1.YSize, SMD1.XSize])
%   MaxOffset: Max offset for which the cross correlation is computed
%              between the two images. 
%              (Pixels)(Default = ceil(SubROISize / 4))
%
% OUTPUT:
%   SubPixelOffsets: The pixel offset of SMD2 relative to SMD1.
%                    (NROIsx2 array)
%   SMDROIs: ROIs of the regions corresponding to the pixel offsets.
%            (NROIsx4 array)([YStart, XStart, YEnd, XEnd])
%   SMDStats: Structure containing some stats about the sub-ROIs.
%
% REQUIRES: 
%   matlab-instrument-control, to use findStackOffset()
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
DataSize = [SMD1.YSize, SMD1.XSize];
if (~exist('SubROISize', 'var') || isempty(SubROISize))
    SubROISize = DataSize.';
end

% Split the images up into the sub-ROIs.
[SubSMDs1, SMDROIs] = smi_helpers.subdivideSMD(SMD1, SubROISize);
SubSMDs2 = smi_helpers.subdivideSMD(SMD2, SubROISize);

% Loop through each ROI and compute the local shift.
NROIs = size(SMDROIs, 1);
SubPixelOffsets = NaN(NROIs, 2);
for nn = 1:NROIs
    % If both SMDs have localizations, estimate the shift.
    if ((numel(SubSMDs1(nn).X)>0) && (numel(SubSMDs1(nn).Y)>0))
        SubPixelOffsets(nn, :) = ...
            smi_core.DriftCorrection.regViaDC(SubSMDs1(nn), SubSMDs2(nn));
    end
end

% If needed, compute some imaging stats. which might be useful to return.
if (nargin > 2)
    SMDStats.NLoc1 = cellfun(@(X) numel(X), {SubSMDs1.X}).';
    SMDStats.NLoc2 = cellfun(@(X) numel(X), {SubSMDs2.X}).';
end


end