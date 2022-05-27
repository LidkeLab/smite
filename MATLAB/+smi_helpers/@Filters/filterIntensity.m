function SMD = filterIntensity(SMD, MeanMultiplier)
%filterIntensity: Filter localizations based on intensity.
%
% INPUTS:
%    SMD              Single Molecule Data structure
%    MeanMultiplier   mean multiplier for which localizations satisfying
%                     intensity > MeanMultiplier*mean(intensity) are removed
%                     [DEFAULT = 2]
%
% OUTPUT:
%    SMD              modified Single Molecule Data structure

% Created by
%    David J. Schodt and Michael J. Wester (5/24/2022)

if ~exist('MeanMultiplier', 'var')
   MeanMultiplier = 2;
end

% Localizations for which intensity > MeanMultiplier*mean(intensity) are
% removed.
MaxPhotons = MeanMultiplier * mean(SMD.Photons);
SMD = ...
   smi_core.SingleMoleculeData.isolateSubSMD(SMD, SMD.Photons <= MaxPhotons);

end
