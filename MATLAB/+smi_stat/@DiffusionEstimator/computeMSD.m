function [MSDEnsemble, MSDIndividual] = computeMSD(TR, MaxFrameLag)
%computeMSD computes the mean squared displacement from TR.
% This method computes the trajectory-wise and ensemble mean squared
% displacements of the trajectories given in 'TR'.
%
% INPUTS:
%   TR: Tracking results structure.
%   MaxFrameLag: Maximum frame difference between localizations to be used
%                in the MSD calculation. (Default = ceil(MaxFrameDiff/4))
%
% OUTPUTS:
%   MSDEnsemble: A structure array containing the ensemble MSD results.
%   MSDIndividual: A structure array containing the trajectory-wise MSD
%                  results.

% Created by:
%   Hanieh Mazloom-Farsibaf (Lidke lab, 2018) in sma-core-alpha
%   rewritten by David J. Schodt (Lidke lab, 2021) in smite


% Set defaults/validate parameters if needed.
MaxFrameDiff = ...
    cell2mat(cellfun(@(X) X(end) - X(1), {TR.FrameNum}, ...
    'UniformOutput', false).');
DefaultMaxFrameLag = ceil(max(MaxFrameDiff) / 4);
if (~exist('MaxFrameLag', 'var') || isempty(MaxFrameLag))
    MaxFrameLag = DefaultMaxFrameLag;
elseif (MaxFrameLag > DefaultMaxFrameLag)
    MaxFrameLag = DefaultMaxFrameLag;
    warning('Input MaxFrameLag=%i is too large.  Default set to %i', ...
        DefaultMaxFrameLag)
end

% Loop through all trajectories in TR and compute the trajectory-wise MSDs.
NTraj = numel(TR);
for ii = 1:NTraj
    
end



















end







