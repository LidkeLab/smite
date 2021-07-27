function [LogLikelihood] = brownianJumpLikelihood(MotionParams, ...
    SquaredDisplacement, FrameLagsAll, LocVarianceSum)
%brownianJumpLikelihood computes the likelihood of given Brownian jumps.
% This method will generate a model of the distribution of the jumps made 
% by a Brownian random walker.
%
% INPUTS:
%   MotionParams: Array of parameters needed for the model. For the one
%                 component model, this is just D. For the two component
%                 model, this is [D_1; D_2; alpha]. For N-components,
%                 this is [D_1, D_2, ..., D_N, 
%                         alpha_1, alpha_2, ..., alpha_{N-1}
%   SquaredDisplacement: The squared displacements made by the
%                        trajectory(ies). (NDatax1 numeric array)
%   FrameLagsAll: All of the frame lags associated with the jumps in
%                 'SquaredDisplacement'. (NDatax1 array)
%   LocVarianceSum: Sum of the localization variances for the two points
%                   used to compute the jumps. This array should be
%                   averaged over x and y. (NJumpsx1 numeric array)
%                   NOTE: I don't know which is better: average the
%                         variances, or average the SEs and square them? My
%                         bet is on averaging variances, but I'm not sure.
%                         This can make a big difference in some cases!
%
% OUTPUTS:
%   LogLikelihood: The log-likelihood of observing the jumps 
%                  'SortedSquaredDisp' assuming they were made by a 
%                   Brownian random walker.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Ensure MotionParams is a column vector.
if isrow(MotionParams)
    MotionParams = MotionParams.';
end

% Determine which model we're using (1-component, 2-component, ...).
NParams = numel(MotionParams);
NComponents = (NParams+1) / 2;

% Pad the MotionParams with the constrained term (i.e., the population
% parameters have to sum to 1 so one of them is defined in terms of the
% others).
MotionParams = [MotionParams; ...
    1 - sum(MotionParams(NComponents+1:NParams))].';

% Compute the likelihood.
% NOTE: I'm taking the mean of the localization variance sums.  That's not
%       ideal, but dealing with those properly becomes too messy/slow
%       (i.e., we get another integral...).
% NOTE: For emphasis, 'LocVarianceSum' is the sum of the localization
%       variances of the two localizations, thus we don't need to
%       multiply it by a factor of 2.
Variance = 2*FrameLagsAll*MotionParams(1:NComponents) ...
    + mean(LocVarianceSum);
LogLikelihood = sum(log(sum((sqrt(SquaredDisplacement)./Variance) ...
    .* exp(-0.5*SquaredDisplacement./Variance) ...
    .* MotionParams((NComponents+1):end), 2)));


end