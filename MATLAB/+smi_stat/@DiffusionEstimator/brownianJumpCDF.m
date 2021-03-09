function [CDFOfJumps] = brownianJumpCDF(...
    MotionParams, SortedJumps, FrameLags, NPoints)
%brownianJumpCDF generates a model of the CDF of Brownian jumps.
% This method will generate a model of the CDF (CPD) of the jumps made by a
% Brownian random walker.
%
% INPUTS:
%   MotionParams: Array of parameters needed for the model. As of now, this
%                 is the variance of the jump sizes for a single frame step
%                 (2 * diffusion constant)
%   SortedJumps: The sorted jumps values used to compute 'CDFOfJumps' in 
%                ascending order. (numeric array)
%   FrameLags: All of the unique frame lags associated with the jumps in
%              'SortedJumps'. Note that this isn't necessarily the same
%              size as SortedJumps, since SortedJumps can contain multiple
%              jumps for each frame lag.
%              (NFrameLagsx1 array)
%   NPoints: The number of data points (or jumps) corresponding to each 
%            frame lag in 'FrameLags' (NFrameLagsx1 array)
%
% OUTPUTS:
%   CDFOfJumps: The CDF (cumulative distribution function, cumulative
%               probability distribution, etc.) of the jump sizes made by a
%               Brownian random walker.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Compute the probability of each frame lag (i.e., the proportion of data
% corresponding to each frame lag).
FrameLagProb = NPoints / sum(NPoints);

% Compute the CDF model.
NFrameLags = numel(FrameLags);
PDFOfJumps = zeros(numel(SortedJumps), 1);
Variance = FrameLags * MotionParams;
for ff = 1:NFrameLags
    PDFOfJumps = PDFOfJumps ...
        + (FrameLagProb(ff) * (SortedJumps/Variance(ff)) ...
        .* exp(-0.5*(SortedJumps.^2)/Variance(ff)));
end
CDFOfJumps = cumsum(PDFOfJumps .* [0; diff(SortedJumps)]);


end