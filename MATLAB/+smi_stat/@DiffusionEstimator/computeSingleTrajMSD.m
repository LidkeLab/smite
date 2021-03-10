function [MSDSingleTraj] = computeSingleTrajMSD(TR, MaxFrameLag, Verbose)
%computeSingleTrajMSD computes the mean squared displacement from TR.
% This method computes the mean squared displacement between localizations
% in the single trajectory provided in TR.
%
% NOTE: This method isn't very memory efficient, as I'm creating two pretty
%       big arrays of zeros (see 'SquaredDisplacement' below) which are
%       usually very sparse. I tried using the sparse matrix built-ins in
%       MATLAB but that made this method too slow.
%
% INPUTS:
%   TR: Tracking results structure containing only one trajectory.  To be
%       sure, only TR(1) will be used in the analysis.
%   MaxFrameLag: Maximum frame difference between localizations to be used
%                in the MSD calculation. (Default is max. possible frame 
%                lag for this trajectory)
%   Verbose: Verbosity level specifying how many temporary outputs should
%            be displayed (e.g., Command Window updates).
%
% OUTPUTS:
%   MSDSingleTraj: A structure array with the following fields:
%       TrajectoryID: The trajectory ID from TR(1).TrajectoryID.
%       MSD: The mean squared displacement between localizations in TR(1)
%            (pixel^2)
%       NCount: The number of displacements used in computing MSD at each
%               point (i.e., MSD(ii) was the mean of NCount(ii) squared
%               displacements).
%       FrameLags: Number of frames between displacements used to compute 
%                  each MSD point (i.e., MSD(ii) is the mean squared
%                  displacement for localizations separated by FrameLag(ii)
%                  frames).
%       SquaredDisplacement: The squared displacements used to compute MSD,
%                            returned as a convenience since it was already
%                            computed internally (i.e., we might as well 
%                            output this if requested).  Each row
%                            corresponds to a single localization, and each
%                            column to a frame lag. (pixel^2)

% Created by:
%   David J. Schodt (Lidke lab, 2021) 
%       based on msdAnalysis.m by Hanieh Mazloom-Farsibaf (Lidke lab, 2018)


% Set defaults if needed.
if (~exist('Verbose', 'var') || isempty(Verbose))
    Verbose = 0;
end
FrameNum = TR(1).FrameNum;
NLocalizations = numel(FrameNum);
DefaultMaxFrameLag = FrameNum(NLocalizations) - FrameNum(1);
if (~exist('MaxFrameLag', 'var') || isempty(MaxFrameLag))
    MaxFrameLag = DefaultMaxFrameLag;
elseif (MaxFrameLag > DefaultMaxFrameLag)
    if (Verbose > 2)
        % For the highest verbosity levels, we should share more info. than
        % just the warning.
        fprintf(['computeSingleTrajMSD(): Input MaxFrameLag=%i but\n', ...
            '\tthe maximum possible frame lag is %i frames.\n', ...
            '\tMaxFrameLag will be set to a default value of\n', ...
            '\tMaxFrameLag=%i\n'], ...
            MaxFrameLag, DefaultMaxFrameLag, DefaultMaxFrameLag)
    end
    MaxFrameLag = DefaultMaxFrameLag;
end

% Loop through localizations and compute the displacement to later
% localizations.
if (Verbose > 2)
    fprintf(['computeSingleTrajMSD(): computing single trajectory\n', ...
        'displacements between %i localizations\n'], NLocalizations)
end
Coordinates = [TR(1).X, TR(1).Y];
MeanVariance = mean([TR(1).X_SE.^2, TR(1).Y_SE.^2], 2);
SquaredDisplacement = zeros(NLocalizations-1, MaxFrameLag);
SquaredDisplacementMask = logical(SquaredDisplacement);
LocVarianceSum = SquaredDisplacement;
for ff = 1:MaxFrameLag
    % Start with the first localization and find the distance to the
    % localization ff frames away.
    Index1 = 1;
    Index2 = 2;
    while (Index2 <= NLocalizations)
        % Check if the current frame gap is that desired.
        FrameDiff = FrameNum(Index2) - FrameNum(Index1);
        if (FrameDiff > ff)
            % This frame gap is too large, so we need a new Index1.
            Index1 = Index1 + 1;
            Index2 = Index1 + 1;
        elseif (FrameDiff == ff)
            % Compute the displacement.
            SquaredDisplacement(Index1, ff) = ...
                (Coordinates(Index1, 1)-Coordinates(Index2, 1))^2 ...
                + (Coordinates(Index1, 2)-Coordinates(Index2, 2))^2;
            SquaredDisplacementMask(Index1, ff) = 1;
            
            % Store the sum of the localization variances.
            LocVarianceSum(Index1, ff) = ...
                MeanVariance(Index1) + MeanVariance(Index2);
            
            % Update the indices, ensuring the gaps don't overlap.
            Index1 = Index2;
            Index2 = Index1 + 1;
        else
            % This frame gap isn't the one desired.
            Index2 = Index2 + 1;
        end
    end
end

% Compute the MSD.
NPoints = sum(logical(SquaredDisplacement), 1).';
KeepBool = logical(NPoints);
if (Verbose > 2)
    fprintf(['computeSingleTrajMSD(): computing single trajectory\n', ...
        'MSD, %i valid frame gaps found\n'], sum(KeepBool))
end
MSD = sum(SquaredDisplacement, 1).' ./ NPoints;
NPoints = NPoints(KeepBool);
MSD = MSD(KeepBool);
FrameLags = (1:MaxFrameLag).';
FrameLags = FrameLags(KeepBool);
MSDSingleTraj.TrajectoryID = TR(1).TrajectoryID;
MSDSingleTraj.MSD = MSD;
MSDSingleTraj.FrameLags = FrameLags;
MSDSingleTraj.NPoints = NPoints;
MSDSingleTraj.SquaredDisplacement = ...
    SquaredDisplacement(SquaredDisplacementMask);
MSDSingleTraj.LocVarianceSum = LocVarianceSum(SquaredDisplacementMask);
FrameLagsAll = FrameLags.' .* SquaredDisplacementMask;
MSDSingleTraj.FrameLagsAll = FrameLagsAll(SquaredDisplacementMask);


end