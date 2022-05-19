function MinMax = setMinMax(SMF)
%setMinMax creates the MinMax structure used throughout the Threshold class
% from a provided SMF structure.  The MinMax fields contain a minimum and a
% maximum value on which to threshold.
%
% INPUT:
%    SMF      Single Molecule Fitting structure
%
% OUTPUT:
%    MinMax   structure containing fields on which to threshold

% Created by
%    Michael J. Wester (2021, Lidke lab)

   MinMax = [];
   MinMax.Photons  = [SMF.Thresholding.MinPhotons, inf];
   MinMax.Bg       = [0, SMF.Thresholding.MaxBg];
   MinMax.PSFSigma = [SMF.Thresholding.MinPSFSigma, ...
       SMF.Thresholding.MaxPSFSigma];
   MinMax.X_SE     = [0, SMF.Thresholding.MaxXY_SE];
   MinMax.Y_SE     = [0, SMF.Thresholding.MaxXY_SE];
   MinMax.Z_SE     = [0, SMF.Thresholding.MaxZ_SE];
   MinMax.PValue   = [SMF.Thresholding.MinPValue, 1];
   MinMax.LogLikelihood = [-inf, inf];

end
