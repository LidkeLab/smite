function SMD = filterIntensity(SMD, Verbose, MeanMultiplier)
%filterIntensity: Filter localizations based on intensity.
%
% INPUTS:
%    SMD              Single Molecule Data structure
%    Verbose          verbosity flag [DEFAULT = false]
%    MeanMultiplier   mean multiplier for which localizations satisfying
%                     intensity > MeanMultiplier*mean(intensity) are removed
%                     [DEFAULT = Inf]
%
% OUTPUT:
%    SMD              modified Single Molecule Data structure

% Created by
%    David J. Schodt and Michael J. Wester (5/24/2022)

if ~exist('Verbose', 'var')
   Verbose = false;
end

if ~exist('MeanMultiplier', 'var')
   MeanMultiplier = Inf;
end

% Localizations for which intensity > MeanMultiplier*mean(intensity) are
% removed.
if ~isinf(MeanMultiplier)
   n_prefilter = numel(SMD.X);
   MaxPhotons = MeanMultiplier * mean(SMD.Photons);
   SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, ...
                                                   SMD.Photons <= MaxPhotons);

   if Verbose >= 2
      fprintf('Intensity filtered localizations kept = %d out of %d\n', ...
              numel(SMD.X), n_prefilter);
   end

end
