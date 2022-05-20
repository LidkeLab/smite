function [SMD] = genTrajGC(SMD, SMF, RhoOff, ...
    NonLinkMarker, UseSparseMatrices)
%genTrajFF connects trajectory segments into longer trajectories.
% This method solves the linear assignment problem for connecting
% trajectory segments into longer trajectories.
%
% INPUTS:
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)
%        containing the localizations that we wish to stitch into
%        trajectories.
%   SMF: Single Molecule Fitting structure defining many of the parameters
%        we'll just to populate cost matrix elements.
%       (see smi_core.SingleMoleculeFitting)
%   RhoOff: Density of dark emitters, given as an image with the same
%           aspect ratio as the SMD coordinate system.
%   NonLinkMarker: A marker in the output CostMatrix that indicates we
%                  don't want to select that element in the linear
%                  assignment problem.
%                  (scalar, ~=inf, ~=nan, and typically < 0)(Default = -1)
%   UseSparseMatrices: This boolean flag will determine whether or not we
%                      should define the gap closing CM as a sparse matrix
%                       (and represent it in MATLAB as a sparse type) or
%                       as a "regular" matrix with a NonLinkMarker for
%                       uninteresting costs.
%                       (boolean flag, 0 or 1)(Default = true)
%
% OUTPUTS:
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)
%        with field ConnectID representing trajectory membership.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults.
if (~exist('NonLinkMarker', 'var') || isempty(NonLinkMarker))
    NonLinkMarker = -1;
end
if (~exist('UseSparseMatrices', 'var') || isempty(UseSparseMatrices))
    UseSparseMatrices = true;
end

% Cluster trajectories together so that we can solve several smaller LAPs.
TrajClusters = smi_cluster.clusterSTDist(SMD, ...
    SMF.Tracking.MaxFrameGap, SMF.Tracking.MaxDistGC);

% Solve the gap-closing LAP(s).
UniqueClusters = unique(TrajClusters);
SMDOut = smi_core.SingleMoleculeData.createSMD();
for nn = UniqueClusters.'
    SMDSub = smi_core.SingleMoleculeData.isolateSubSMD(SMD, ...
        TrajClusters == nn);
    SMDSub.ConnectID = smi_helpers.compressToRange(SMDSub.ConnectID);
    CostMatrix = smi.SPT.createCostMatrixGC(SMDSub, SMF, RhoOff, ...
        NonLinkMarker, UseSparseMatrices);
    Link12 = smi.SPT.solveLAP(CostMatrix);
    SMDSub = smi.SPT.connectTrajGC(SMDSub, Link12);
    SMDOut = smi_core.SingleMoleculeData.catSMD(SMDOut, SMDSub, ...
        false, false);
end
SMD = SMDOut;


end