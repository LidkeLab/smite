function [MSDSingleTraj] = computeSingleTrajMSD(TR, FrameLagRange, Verbose)
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
%   FrameLagRange: Range of frame differences between localizations for
%                  which MSD is computed. (Default is from 1 to 1/4 of the 
%                  max possible frame lag for this trajectory)
%                  ([min., max.])
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
if (~exist('FrameLagRange', 'var') || isempty(FrameLagRange))
    FrameLagRange = [1, FrameNum(NLocalizations) - FrameNum(1)];
end

% Loop through localizations and compute the displacement to later
% localizations.
if (Verbose > 2)
    fprintf(['computeSingleTrajMSD(): computing single trajectory\n', ...
        'displacements between %i localizations\n'], NLocalizations)
end
Coordinates = [TR(1).X, TR(1).Y];
MeanVariance = mean([TR(1).X_SE.^2, TR(1).Y_SE.^2], 2);
FrameLags = (FrameLagRange(1):FrameLagRange(2)).';
NFrameLags = numel(FrameLags);
SquaredDisplacement = zeros(NLocalizations-1, NFrameLags);
SquaredDisplacementMask = logical(SquaredDisplacement);
LocVarianceSum = SquaredDisplacement;
for ff = 1:NFrameLags
    % Start with the first localization and find the distance to the
    % localization ff frames away.
    Index1 = 1;
    Index2 = 2;
    while (Index2 <= NLocalizations)
        % Check if the current frame gap is that desired.
        FrameDiff = FrameNum(Index2) - FrameNum(Index1);
        if (FrameDiff > FrameLags(ff))
            % This frame gap is too large, so we need a new Index1.
            Index1 = Index1 + 1;
            Index2 = Index1 + 1;
        elseif (FrameDiff == FrameLags(ff))
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
FrameLagsAll = FrameLags.' .* SquaredDisplacementMask;
FrameLags = FrameLags(KeepBool);
MSDSingleTraj.ConnectID = TR(1).ConnectID;
MSDSingleTraj.MSD = MSD;
MSDSingleTraj.FrameLags = FrameLags;
MSDSingleTraj.FrameLagsAll = FrameLagsAll(SquaredDisplacementMask);
MSDSingleTraj.NPoints = NPoints;
MSDSingleTraj.SquaredDisplacement = ...
    SquaredDisplacement(SquaredDisplacementMask);
MSDSingleTraj.LocVarianceSum = LocVarianceSum(SquaredDisplacementMask);


end