function [TransitionMatrix] = ...
    generateTransitionMatrix(TransitionRates, DeltaT)
%generateTransitionMatrix generates the transition matrix of a Markov model
% This method will generate a transition matrix for a Markov model given
% the input TransitionRates.

% INPUTS:
%   TransitionRates: Array containing the transition rate constants
%                    between the states of the Markov model.  The element
%                    TransitionRates(i, j) is the transition rate from
%                    state i to state j.  TransitionRates(i, i) is not
%                    used for any value of i.
%                    (NxN, N = number of states)(1 / frames)
%   DeltaT: Time step(s) between observations.
%           (1x1 or NObservations-1 x 1 array)(frames)
%
% OUTPUTS:
%   TransitionMatrix: Matrix containing the transition probabilities
%                     between the N states of the Markov model.
%                     TransitionMatrix(i, j, k) gives the probability of
%                     transition from state i to state j during DeltaT(k).
%                     Note that sum(TransitionMatrix, 2, ii) = ones(N, 1)
%                     (NStates x NStates x numel(DeltaT))

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Force the inputs to be double precision floats (just in case).
TransitionRates = double(TransitionRates);
DeltaT = double(DeltaT);

% Ensure diagonal terms of the input TransitionRates matrix are 0.
NStates = size(TransitionRates, 1);
NObservations = numel(DeltaT);
DiagonalBool = logical(eye(NStates));
TransitionRates(DiagonalBool) = zeros(NStates, 1);

% Loop through each element in DeltaT and compute the corresponding
% transition matrix.
LeavingProbability = repmat(sum(TransitionRates, 2), ...
    [1, NStates]);
TransitionMatrix = zeros(NStates, NStates, NObservations);
for ii = 1:NObservations
    % Populate the off-diagonal terms of the transition matrix.  These
    % terms correspond to transitions to other states.
    ArrivalProbability = 1 - exp(-LeavingProbability*DeltaT(ii));
    TransitionMatrixCurrent = TransitionRates .* ArrivalProbability ...
        ./ LeavingProbability;
    
    % Populate the diagonal elements of the transition matrix.  These terms
    % correspond to remaining in the same state.
    % NOTE: The diagonal elements are computed to satisfy
    %       sum(TransitionMatrix, 2) = ones(NStates, 1)
    TransitionMatrixCurrent(DiagonalBool) = ...
        1 - sum(TransitionMatrixCurrent, 2);
    
    % Add the current TransitionMatrix to the full TransitionMatrix array.
    TransitionMatrix(:, :, ii) = TransitionMatrixCurrent;
end


end