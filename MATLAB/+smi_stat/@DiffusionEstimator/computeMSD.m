function [MSDSingleTraj, MSDEnsemble] = ...
    computeMSD(TR, FrameLagRange, Verbose)
%computeMSD computes the mean squared displacement from TR.
% This method computes the trajectory-wise and ensemble mean squared
% displacements of the trajectories given in 'TR'.
%
% INPUTS:
%   TR: Tracking results structure.
%   FrameLagRange: Range of frame differences between localizations for
%                  which MSD is computed. (Default is from 1 to 1/4 of the 
%                  max possible frame lag for the ensemble
%                  calculation)([min., max.])
%   Verbose: Verbosity level specifying how many temporary outputs should
%            be displayed (e.g., Command Window updates).
%
% OUTPUTS:
%   MSDSingleTraj: A structure array containing the trajectory-wise MSD
%                  results.  All units are in camera units (pixels, frames)
%   MSDEnsemble: A structure array containing the ensemble MSD results.
%                All units are in camera units (pixels, frames)

% Created by:
%   David J. Schodt (Lidke lab, 2021)
%       based on msdAnalysis.m by Hanieh Mazloom-Farsibaf (Lidke lab, 2018)


% Set defaults/validate parameters if needed.
if (~exist('Verbose', 'var') || isempty(Verbose))
    Verbose = 0;
end
MaxFrameDiff = ...
    max(cell2mat(cellfun(@(X) X(end) - X(1), {TR.FrameNum}, ...
    'UniformOutput', false).'));
DefaultMaxFrameLag = ceil(MaxFrameDiff / 4);
if (~exist('FrameLagRange', 'var') || isempty(FrameLagRange))
    FrameLagRange = [1, DefaultMaxFrameLag];
elseif (FrameLagRange(2) > MaxFrameDiff)
    if (Verbose > 2)
        % For the highest verbosity levels, we should share more info. than
        % just the warning.
        fprintf(['computeMSD(): Input max(FrameLagRange)=%i but the\n', ...
            '\tmaximum possible frame lag is %i frames. FrameLagRange\n', ...
            '\twill be set to a default value of [1, %i]\n'], ...
            FrameLagRange(2), MaxFrameDiff, DefaultMaxFrameLag)
    elseif (Verbose > 0)
        warning(['computeMSD(): Input max(FrameLagRange)=%i is too large. ', ...
            'Using default of %i.'], ...
            FrameLagRange(2), DefaultMaxFrameLag)
    end
    FrameLagRange = [1, DefaultMaxFrameLag];
end

% Loop through all trajectories in TR and compute the trajectory-wise and
% ensemble MSDs.
NTraj = numel(TR);
if (Verbose > 1)
    fprintf(['computeMSD(): computing trajectory-wise MSDs for %i ', ...
        'localizations.\n'], NTraj)
end
MSDSingleTraj = struct([]);
MSDMatrix = zeros(NTraj, FrameLagRange(2));
NPointsMatrix = MSDMatrix;
SquaredDisplacement = [];
LocVarianceSum = [];
FrameLagsAll = [];
for ii = 1:NTraj
    % Compute the MSD for this trajectory.
    if (Verbose > 2)
        fprintf(['computeMSD(): computing MSD for ', ...
            'trajectory TR(%i)...\n'], ii)
    end
    MSDCurrent = smi_stat.DiffusionEstimator.computeSingleTrajMSD(...
        TR(ii), FrameLagRange, Verbose);
    MSDSingleTraj = [MSDSingleTraj; MSDCurrent];
    
    % Store the single trajectory MSD in a matrix with all of the
    % trajectory MSDs.
    CurrentLags = MSDCurrent.FrameLags;
    MSDMatrix(ii, CurrentLags) = MSDCurrent.MSD;
    NPointsMatrix(ii, CurrentLags) = MSDCurrent.NPoints;
    SquaredDisplacement = [SquaredDisplacement; ...
        MSDCurrent.SquaredDisplacement];
    LocVarianceSum = [LocVarianceSum; MSDCurrent.LocVarianceSum];
    FrameLagsAll = [FrameLagsAll; MSDCurrent.FrameLagsAll];
end
if (Verbose > 1)
    fprintf('computeMSD(): computing ensemble MSD...\n')
end
FrameLags = (FrameLagRange(1):FrameLagRange(2)).';
MSDMatrix = MSDMatrix(:, FrameLags);
NPointsMatrix = NPointsMatrix(:, FrameLags);
NPoints = sum(NPointsMatrix, 1).';
MSD = sum(MSDMatrix.*NPointsMatrix, 1).' ./ NPoints;
KeepBool = ~isnan(MSD);
MSDEnsemble.MSD = MSD(KeepBool);
MSDEnsemble.FrameLags = FrameLags(KeepBool);
MSDEnsemble.NPoints = NPoints(KeepBool);
MSDEnsemble.SquaredDisplacement = SquaredDisplacement;
MSDEnsemble.LocVarianceSum = LocVarianceSum;
MSDEnsemble.FrameLagsAll = FrameLagsAll;


end