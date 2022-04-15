function [SMD, SMDMasked, KeepBool] = maskSMD(SMD, Mask)
%maskSMD masks the input SMD based on the image Mask.
% This method will isolate a subset of SMD defined by the boolean mask
% defined as an image in 'Mask'.
%
% INPUTS:
%   SMD: A Single Molecule Data structure.
%   Mask: An SMD.YSize x SMD.XSize boolean image, where false values
%         represent regions that should be masked in the SMD.  If the mask
%         image is a different size than needed for SMD, it will be
%         silently resized within this method!
%
% OUTPUTS:
%   SMD: A Single Molecule Data structure containing localizations from SMD 
%        that fall within true pixels of 'Mask'.
%   SMDMasked: Localizations of the input 'SMD' that fell within the false
%              pixels of 'Mask'.
%   KeepBool: Boolean array used to define the output 'SMD' from the input
%             'SMD'.

% Created by:
%   David J. Schodt (Lidke lab, 2022)


% Define a boolean array for masking the SMD, accounting for the mask
% magnification if needed (e.g., we might have an upsampled mask to get
% sub-pixel masking resolution).
KeepBool = smi_core.SingleMoleculeData.defineSMDMask(SMD, Mask);

% Extract the desired SMD as well as the masked localizations.
SMDMasked = smi_core.SingleMoleculeData.isolateSubSMD(SMD, ~KeepBool);
SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, KeepBool);


end