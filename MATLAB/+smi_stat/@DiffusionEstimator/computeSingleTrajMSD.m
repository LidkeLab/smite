function [MSDSingleTraj] = computeSingleTrajMSD(TR)
%computeSingleTrajMSD computes the mean squared displacement from TR.
% This method computes the mean squared displacement between localizations
% in the single trajectory provided in TR.
%
% NOTE: Enforcing a max. frame lag is best done outside of this method.  In
%       testing, enforcing that max. lag in this code showed negligible
%       speed improvements, thus it didn't seem beneficial to throw out the
%       large lag data internally.
%
% INPUTS:
%   TR: Tracking results structure containing only one trajectory.  To be
%       sure, only TR(1) will be used in the analysis.
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


% Loop through localizations and compute the displacement to later
% localizations.
Coordinates = [TR(1).X, TR(1).Y];
FrameNum = TR(1).FrameNum;
NLocalizations = numel(FrameNum);
MaxFrameLag = FrameNum(NLocalizations) - FrameNum(1);
SquaredDisplacement = zeros(NLocalizations-1, MaxFrameLag);
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
FrameLags = (1:MaxFrameLag).';
NPoints = sum(logical(SquaredDisplacement), 1).';
MSD = sum(SquaredDisplacement, 1).' ./ NPoints;
KeepBool = ~isnan(MSD);
NPoints = NPoints(KeepBool);
MSD = MSD(KeepBool);
FrameLags = FrameLags(KeepBool);
MSDSingleTraj.TrajectoryID = TR(1).TrajectoryID;
MSDSingleTraj.MSD = MSD;
MSDSingleTraj.FrameLags = FrameLags;
MSDSingleTraj.NPoints = NPoints;
MSDSingleTraj.SquaredDisplacement = SquaredDisplacement;


end