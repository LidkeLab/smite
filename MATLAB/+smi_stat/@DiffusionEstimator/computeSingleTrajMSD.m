function [MSD, NCount, SquaredDisplacement] = computeSingleTrajMSD(TR)
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
%   MSD: The mean squared displacement between localizations in TR(1).
%   NCount: The number of displacements used in computing MSD at each point
%           (i.e., MSD(ii) was the mean of NCount(ii) squared
%           displacements).
%   SquaredDisplacement: The squared displacements used to compute MSD,
%                        returned as a convenience since it was already
%                        computed internally (i.e., we might as well output
%                        this if requested).

% Created by:
%   Hanieh Mazloom-Farsibaf (Lidke lab, 2018) in msdAnalysis.m
%   rewritten by David J. Schodt (Lidke lab, 2021) in smite


% Loop through localizations and compute the displacement to later
% localizations.
Coordinates = [TR(1).X, TR(1).Y];
FrameNum = TR(1).FrameNum;
NPoints = numel(FrameNum);
SquaredDisplacement = zeros(NPoints * (NPoints-1) / 2, 1);
FrameDiff = SquaredDisplacement;
for ii = 1:(NPoints-1)
    for jj = ii+1:NPoints
        SquaredDisplacement(ii + jj*(ii-1)) = ...
            (Coordinates(ii, 1)-Coordinates(jj, 1))^2 ...
            + (Coordinates(ii, 2)-Coordinates(jj, 2))^2;
        FrameDiff(ii + jj*(ii-1)) = FrameNum(jj) - FrameNum(ii);     
    end
end

% Compute the MSD.
MaxFrameLag = max(FrameDiff);
MSD = NaN(MaxFrameLag, 1);
NPoints = MSD;
for ff = 1:MaxFrameLag
    CurrentLagBool = (FrameDiff == ff);
    NPoints(ff) = sum(CurrentLagBool);
    MSD(ff) = sum(SquaredDisplacement(CurrentLagBool)) / NPoints(ff);
end


end