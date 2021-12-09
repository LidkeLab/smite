function [SMD] = genTrajFF(SMD, SMF, RhoOff, NonLinkMarker)
%genTrajFF connects localizations frame-to-frame into trajectories.
% This method loops through frames of localizations in SMD and stitches
% localizations into trajectories frame-to-frame.
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

% Loop over frames and solve the frame-to-frame connection LAP.
SMD.ConnectID = (1:numel(SMD.FrameNum)).';
UniqueFrames = unique(SMD.FrameNum, 'sorted');
for ff = UniqueFrames(1:(end-1)).'
    % Create the frame-to-frame connection cost matrix.
    CostMatrix = smi.SPT.createCostMatrixFF(SMD, SMF, RhoOff, ...
        ff, NonLinkMarker);
    if (numel(CostMatrix) < 2)
        % If there's only one localization considered, there's no use in
        % proceeding.
        continue
    end
    
    % Perform the linear assignment problem to determine how we should link
    % together trajectories.
    Link12 = smi.SPT.solveLAP(CostMatrix);
    SMD = smi.SPT.connectTrajFF(SMD, Link12, ff);
end


end