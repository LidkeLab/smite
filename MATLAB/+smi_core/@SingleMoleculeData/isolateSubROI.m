function [SMD] = isolateSubROI(SMD, ROI)
%isolateSubROI isolates localizations within ROI.
% This method will isolate a subset of SMD defined by the region of
% interest ROI.
%
% NOTE: As written, the pixel convention used in this method is that the
%       left edge of pixel n has coordinate n-1 and the right edge has
%       coordinate n (where entries of ROI are defined as integers n>0).
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


KeepBool = ((SMD.Y>=(ROI(1)-1)) & (SMD.Y<=ROI(3)) ...
    & (SMD.X>=(ROI(2)-1)) & (SMD.X<=ROI(4)));
SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, KeepBool);


end