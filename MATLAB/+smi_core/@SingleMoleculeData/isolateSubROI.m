function [SMD] = isolateSubROI(SMD, ROI)
%isolateSubROI isolates localizations within ROI.
% This method will isolate a subset of SMD defined by the region of
% interest ROI.
%
% INPUTS:
%   SMD: A Single Molecule Data structure.
%   ROI: Region of interested with inclusive edges. 
%        (pixels)([YStart, XStart, YEnd, XEnd])
%
% OUTPUTS:
%   SMD: A Single Molecule Data structure whose vector fields correspond to
%        the input region of interested defined by ROI.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


KeepBool = ((SMD.Y>=ROI(1)) & (SMD.Y<=ROI(3)) ...
    & (SMD.X>=ROI(2)) & (SMD.X<=ROI(4)));
SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, KeepBool);


end