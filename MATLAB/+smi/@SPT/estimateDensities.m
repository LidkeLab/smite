function [RhoOff, RhoOn] = estimateDensities(SMD, SMF)
%estimateDensities estimates emitter densities.
% This method will make an estimate of the density of emitters based on
% trajectories defined by SMD (i.e., SMD with field ConnectID populated to
% indicate trajectory membership).
%
% INPUTS:
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)
%        containing localizations associated by SMD.ConnectID.
%   SMF: Single Molecule Fitting structure (see
%        smi_core.SingleMoleculeFitting).
%
% OUTPUTS:
%   RhoOff: Density of dark emitters (emitters capable of transitioning to
%           a fluorescent state). (emitters / pixel^2)
%   RhoOn: Density of visible emitters. (emitters / pixel^2)
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Estimate the average density of visible emitters.
RhoOn = mean(smi_core.SingleMoleculeData.computeDensity(SMD));

% Estimate the density of dark emitters.
RhoOff = RhoOn * (SMF.Tracking.K_off/SMF.Tracking.K_on);


end