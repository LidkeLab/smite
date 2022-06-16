function SMD = inflateSE(SMD, Verbose, SEAdjust)
%inflateSE: Inflate standard errors.
%
% INPUTS:
%    SMD        Single Molecule Data structure
%    Verbose    verbosity flag [Default = false]
%    SEAdjust   standard error inflation amount (pixel) [DEFAULT = 0]
%
% OUTPUT:
%    SMD        modified Single Molecule Data structure

% Created by
%    David J. Schodt and Michael J. Wester (5/24/2022)

if ~exist('Verbose', 'var')
   Verbose = false;
end

if ~exist('SEAdjust', 'var')
   SEAdjust = 0;
end

if SEAdjust > 0
   SMD.X_SE = SMD.X_SE + SEAdjust;
   SMD.Y_SE = SMD.Y_SE + SEAdjust;

   if Verbose >= 2
      fprintf('Inflated standard errors by %g pixels.\n', SEAdjust);
   end
end

end
