function [MSDSingleTraj, MSDEnsemble] = ...
    computeMSD(TR, MaxFrameLag, Verbose)
%computeMSD computes the mean squared displacement from TR.
% This method computes the trajectory-wise and ensemble mean squared
% displacements of the trajectories given in 'TR'.
%
% INPUTS:
%   TR: Tracking results structure.
%   MaxFrameLag: Maximum frame difference between localizations to be used
%                in the MSD calculation. (Default is 1/4 of the max 
%                possible frame lag for the ensemble calculation)
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
if (~exist('MaxFrameLag', 'var') || isempty(MaxFrameLag))
    MaxFrameLag = DefaultMaxFrameLag;
elseif (MaxFrameLag > DefaultMaxFrameLag)
    if (Verbose > 2)
        % For the highest verbosity levels, we should share more info. than
        % just the warning.
        fprintf(['computeMSD(): Input MaxFrameLag=%i but the maximum\n',...
            '\tpossible frame lag is %i frames. MaxFrameLag will be\n', ...
            '\tset to a default value of MaxFrameLag=%i\n'], ...
            MaxFrameLag, MaxFrameDiff, DefaultMaxFrameLag)
    elseif (Verbose > 0)
        warning(['computeMSD(): Input MaxFrameLag = %i is too large. ', ...
            'Using default of %i.'], ...
            MaxFrameLag, DefaultMaxFrameLag)
    end
    MaxFrameLag = DefaultMaxFrameLag;
end

% Loop through all trajectories in TR and compute the trajectory-wise and
% ensemble MSDs.
NTraj = numel(TR);
if (Verbose > 1)
    fprintf(['computeMSD(): computing trajectory-wise MSDs for %i ', ...
        'localizations.\n'], NTraj)
end
MSDSingleTraj = struct([]);
MSDMatrix = zeros(NTraj, MaxFrameLag);
NPointsMatrix = MSDMatrix;
for ii = 1:NTraj
    % Compute the MSD for this trajectory.
    if (Verbose > 2)
        fprintf(['computeMSD(): computing MSD for ', ...
            'trajectory TR(%i)...\n'], ii)
    end
    MSDCurrent = smi_stat.DiffusionEstimator.computeSingleTrajMSD(...
        TR(ii), MaxFrameLag, Verbose);
    MSDSingleTraj = [MSDSingleTraj; MSDCurrent];
    
    % Store the single trajectory MSD in a matrix with all of the
    % trajectory MSDs.
    CurrentLags = MSDCurrent.FrameLags;
    MSDMatrix(ii, CurrentLags) = MSDCurrent.MSD;
    NPointsMatrix(ii, CurrentLags) = MSDCurrent.NPoints;
end
if (Verbose > 1)
    fprintf('computeMSD(): computing ensemble MSD...\n')
end
FrameLags = (1:MaxFrameLag).';
MSDMatrix = MSDMatrix(:, FrameLags);
NPointsMatrix = NPointsMatrix(:, FrameLags);
NPoints = sum(NPointsMatrix, 1).';
MSD = sum(MSDMatrix.*NPointsMatrix, 1).' ./ NPoints;
KeepBool = ~isnan(MSD);
MSDEnsemble.MSD = MSD(KeepBool);
MSDEnsemble.FrameLags = FrameLags(KeepBool);
MSDEnsemble.NPoints = NPoints(KeepBool);


end