function SMD = filterFC(SMD, Verbose, nFC)
%filterFC: Filter out localizations representing fewer than nFC frame connections.
%
% INPUTS:
%    SMD     Single Molecule Data structure
%    Verbose verbosity flag [DEFAULT = false]
%    nFC     minimum number of frame connected localizations representing a
%            single localization allowed [DEFAULT = 1]
%
% OUTPUT:
%    SMD     modified Single Molecule Data structure

% Created by
%    David J. Schodt and Michael J. Wester (5/24/2022)

if ~exist('Verbose', 'var')
   Verbose = false;
end

if ~exist('nFC', 'var')
   nFC = 1;
end

% Remove localizations representing fewer than nFC frame-connected
% localizations.
if nFC > 1
   n_prefilter = numel(SMD.X);
   SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, SMD.NCombined >= nFC);

   if Verbose >= 2
      fprintf(['Frame connected filtered localizations kept = %d ', ...
               'out of %d\n'], numel(SMD.X), n_prefilter);
   end
end

end
