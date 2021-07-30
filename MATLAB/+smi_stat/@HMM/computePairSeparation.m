function [TRArray] = computePairSeparation(TRArray)
%computePairSeparation computes separation between traj. pairs in TRArray.
% This method will loop over the trajectory pairs in TRArray, compute the
% center-to-center distance between those trajectories, and then store the
% result in the output TRArray.  This method also computes the x,y averaged
% standard errors (this is a bit out of place based on this methods name,
% but it seemed to make sense to stick it in here since the x,y averaged
% standard errors are almost always something we'll need when we look at
% the separation between two trajectories).
%
% NOTE: This method assumes that the pairs in TRArray are already matched
%       in frame number, i.e. 
%       all(TRArray(1, ii).FrameNum == TRArray(2, ii).FrameNum) == 1
%
% INPUTS:
%   TRArray: A structure array of TR structures containing dimer candidates 
%            (see HMM.findDimerCandidates()) which have been "matched"
%            using HMM.isolateCandidateTRArray().
%
% OUTPUTS:
%   TRArray: Input TRArray but with additional fields Separation
%            (separation between the two trajectories in a column of
%            TRArray) and AverageSE (the averaged x,y averaged standard
%            error for a given trajectory).
%           

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Loop through the TRArray and update it to contain the Separation and
% AverageSE fields.
for ii = 1:size(TRArray, 2)
    XTraj1 = double(TRArray(1, ii).X);
    XSETraj1 = double(TRArray(1, ii).X_SE);
    YTraj1 = double(TRArray(1, ii).Y);
    YSETraj1 = double(TRArray(1, ii).Y_SE);
    XTraj2 = double(TRArray(2, ii).X);
    XSETraj2 = double(TRArray(2, ii).X_SE);
    YTraj2 = double(TRArray(2, ii).Y);
    YSETraj2 = double(TRArray(2, ii).Y_SE);
    Separation = sqrt((XTraj1-XTraj2).^2 + (YTraj1-YTraj2).^2);
    AverageSETraj1 = mean([XSETraj1, YSETraj1], 2);
    AverageSETraj2 = mean([XSETraj2, YSETraj2], 2);
    TRArray(1, ii).Separation = Separation;
    TRArray(1, ii).AverageSE = AverageSETraj1;
    TRArray(2, ii).Separation = Separation;
    TRArray(2, ii).AverageSE = AverageSETraj2;
end


end