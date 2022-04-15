function [KeepBool] = defineSMDMask(SMD, Mask)
%defineSMDMask defines a boolean for masking SMD localizaions.
% This method takes the boolean image 'Mask' and converts it to a boolean
% array to mask the provided SMD.  That is, the output 'KeepBool' indicates
% localizations in 'SMD' that fall within the hot bits of 'Mask'.
%
% INPUTS:
%   SMD: A Single Molecule Data structure.
%   Mask: An SMD.YSize x SMD.XSize boolean image, where false values
%         represent regions that should be masked in the SMD.  If the mask
%         image is a different size than needed for SMD, it will be
%         silently resized within this method!
%
% OUTPUTS:
%   KeepBool: Boolean array used to define the output 'SMD' from the input
%             'SMD'.

% Created by:
%   David J. Schodt (Lidke lab, 2022)


% Validate the size of the input mask and define a rescaling if needed.
MaskSize = size(Mask);
Mag = max([SMD.YSize, SMD.XSize] ./ MaskSize);

% Define a boolean array for masking the SMD, accounting for the mask
% magnification if needed (e.g., we might have an upsampled mask to get
% sub-pixel masking resolution).
XInds = min(SMD.XSize, max(1, round(Mag*(SMD.X-0.5) + 0.5)));
YInds = min(SMD.YSize, max(1, round(Mag*(SMD.Y-0.5) + 0.5)));
KeepBool = Mask(sub2ind([SMD.YSize, SMD.XSize], YInds, XInds));


end