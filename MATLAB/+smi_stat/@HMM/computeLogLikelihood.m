function [LogLikelihood] = ...
    computeLogLikelihood(TransitionMatrixSeries, EmissionProbabilitySeries)
%computeLogLikelihood computes the LogLikelihood of a sequence in the HMM.
% This method will compute the log likelihood of observing a pair of
% trajectories given a specified Hidden Markov Model (HMM).

% INPUTS:
%   TransitionMatrixSeries: An array containing transition matrices for
%                           each observation in the sequence. 
%                           (NStates x NStates x NObservations-1 double)
%   EmissionProbabilitySeries: An array containing the emission
%                              probabilities for each state at each
%                              observation. 
%                              (NStates x NObservations double)
%
% OUTPUTS:
%   LogLikelihood: The natural logarithm of the likelikhood of observing
%                  the input data given the specified HMM.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Loop through the observations and compute the log-likelihood.
NObservations = size(EmissionProbabilitySeries, 1);
LikelihoodArray = EmissionProbabilitySeries(1, :);
LikelihoodScaling = zeros(NObservations, 1);
LikelihoodScaling(1) = sum(LikelihoodArray);
for ii = 2:NObservations
    % Update the likelihood matrix.
    % NOTE: You might see the emission probabilities forced into a diagonal
    %       matrix in some literature, but we don't need to do that.  It's
    %       in an array form as EmissionProbabilitySeries(ii, :) so I just
    %       use element-wise multiplication to get the same effect with
    %       fewer computations needed.
    LikelihoodArray = LikelihoodArray ...
        * TransitionMatrixSeries(:, :, ii-1) ...
        .* EmissionProbabilitySeries(ii, :);
    LikelihoodScaling(ii) = sum(LikelihoodArray);
    LikelihoodArray = LikelihoodArray ./ LikelihoodScaling(ii);
end
LogLikelihood = sum(log(LikelihoodScaling));


end