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

   SMD = filterNonNeg(SMD)
   SMD = filterIntensity(SMD, MeanMultiplier)
   SMD = inflateSE(SMD, SEAdjust)
   SMD = filterFC(SMD, nFC)
   SMD = filterNN(SMD, n_NN, MedianMultiplier)

end % methods(Static)

end % classdef Filters
