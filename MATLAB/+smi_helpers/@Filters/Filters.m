classdef Filters < handle

% Filters operating on SMD (Single Molecule Data) structures used especially
% for BaGoL analyses:
%
% SR data -> remove localizations with negative coordinates
%         -> intensity filter
%         -> inflate standard errors
%         -> frame connection, removing connections which involve only 1 frame
%         -> NND filter --- Do not use on dSTORM data!
%         -> BaGoL

methods(Static)

   SMD = filterNonNeg(SMD, Verbose)
   SMD = filterIntensity(SMD, Verbose, MeanMultiplier)
   SMD = inflateSE(SMD, Verbose, SEAdjust)
   SMD = filterFC(SMD, Verbose, nFC)
   SMD = filterNN(SMD, Verbose, n_NN, MedianMultiplier)

   SMD = filterImag(SMD, Verbose)

end % methods(Static)

end % classdef Filters
