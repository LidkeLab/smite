function [CDFOfJumps] = brownianJumpCDF(MotionParams, ...
    SortedSquaredDisp, FrameLags, NPoints, LocVarianceSum)
%brownianJumpCDF generates a model of the CDF of Brownian jumps.
% This method will generate a model of the CDF (CPD) of the jumps made by a
% Brownian random walker.
%
% INPUTS:
%   MotionParams: Array of parameters needed for the model. For the one
%                 component model, this is just D. For the two component
%                 model, this is [D_1; D_2; alpha_1, alpha_2].
%                 For N-components, this is [D_1, D_2, ..., D_N, 
%                         alpha_1, alpha_2, ..., alpha_N
%   SortedSquaredDisp: The sorted (in ascending order) squared jumps used 
%                      to compute 'CDFOfJumps'. (numeric array)
%   FrameLags: All of the unique frame lags associated with the jumps in
%              'SortedSquaredDisp'. Note that this isn't necessarily the 
%              same size as 'SortedSquaredDisp', since 'SortedSquaredDisp'
%              can contain multiple jumps for each frame lag. 
%              (NFrameLagsx1 array)
%   NPoints: The number of data points (or jumps) corresponding to each
%            frame lag in 'FrameLags' (NFrameLagsx1 array)
%   LocVarianceSum: Sum of the localization variances for the two points
%                   used to compute the jumps. This array should be
%                   averaged over x and y.
%                   (NJumpsx1 numeric array)
%                   NOTE: I don't know which is better: average the
%                         variances, or average the SEs and square them? My
%                         bet is on averaging variances, but I'm not sure.
%                         This can make a big difference in some cases!
%
% OUTPUTS:
%   CDFOfJumps: The CDF (cumulative distribution function, cumulative
%               probability distribution, etc.) of the jump sizes made by a
%               Brownian random walker.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Ensure MotionParams is a column vector.
if isrow(MotionParams)
    MotionParams = MotionParams.';
end

% Determine which model we're using (1-component, 2-component, ...).
NParams = numel(MotionParams);
NComponents = NParams / 2;

% Compute the probability of each frame lag (i.e., the proportion of data
% corresponding to each frame lag).
FrameLagProb = NPoints / sum(NPoints);

% Compute the CDF model.
% NOTE: I'm taking the mean of the localization variance sums.  That's not
%       ideal, but dealing with those properly becomes too messy/slow
%       (i.e., we get another integral...).
NFrameLags = numel(FrameLags);
ZerosArray = zeros(numel(SortedSquaredDisp), 1);
CDFOfJumps = ZerosArray;
for nn = 1:NComponents
    % Define the variance term for this component.
    % NOTE: For emphasis, 'LocVarianceSum' is the sum of the localization
    %       variances of the two localizations, thus we don't need to
    %       multiply it by a factor of 2.
    Variance = 2*FrameLags*MotionParams(nn) + mean(LocVarianceSum);
    
    % Sum over frame lags.
    CDFOfJumpsNN = ZerosArray;
    for ff = 1:NFrameLags
        CDFOfJumpsNN = CDFOfJumpsNN ...
            + (FrameLagProb(ff) ...
            * (1-exp(-0.5*SortedSquaredDisp/Variance(ff))));
    end
    CDFOfJumps = CDFOfJumps + MotionParams(nn+NComponents)*CDFOfJumpsNN;
end


end