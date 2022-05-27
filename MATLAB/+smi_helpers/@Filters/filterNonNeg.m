function SMD = filterNonNeg(SMD)
%filterNonNeg: Filter out localizations with negative coordinates.
%
% INPUT:
%    SMD   Single Molecule Data structure
%
% OUTPUT:
%    SMD   modified Single Molecule Data structure

% Created by
%    David J. Schodt and Michael J. Wester (5/24/2022)

SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, ...
         SMD.X >= 0 & SMD.Y >= 0);

end
