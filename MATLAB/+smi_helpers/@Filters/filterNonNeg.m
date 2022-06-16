function SMD = filterNonNeg(SMD, Verbose)
%filterNonNeg: Filter out localizations with negative coordinates.
%
% INPUT:
%    SMD     Single Molecule Data structure
%    Verbose verbosity flag [DEFAULT = false]
%
% OUTPUT:
%    SMD   modified Single Molecule Data structure

% Created by
%    David J. Schodt and Michael J. Wester (5/24/2022)

if ~exist('Verbose', 'var')
   Verbose = false;
end

n_prefilter = numel(SMD.X);
SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, SMD.X >= 0 & SMD.Y >= 0);

if Verbose >= 2
   fprintf('Nonnegative localizations kept = %d out of %d\n', ...
           numel(SMD.X), n_prefilter);
end

end
