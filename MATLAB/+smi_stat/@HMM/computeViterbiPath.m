function [StateSequence] = computeViterbiPath(...
    StateSpace, InitialProbability, TransitionMatrix, EmissionMatrix)
%computeViterbiPath estimates state sequence for a HMM using Viterbi alg.
% This method uses the Viterbi algorithm to estimate the most likely state
% sequence of hidden states in a hidden Markov model (HMM) given a set of
% observations.  In this method, the observations are intrinsically 
% accounted for in the EmissionMatrix.
% 
% NOTE: The observation space is assumed to be infinite, e.g., as it would
%       be if the observation space is the set of all possible distances 
%       between two trajectories.  In this case, the columns of the
%       emission matrix MUST correspond directly to an observation given in
%       Observation, i.e., the emission probabilities for Obervation(ii)
%       are given in the column EmissionMatrix(:, ii).
% 
% INPUTS:
%   StateSpace: The set of states in the HMM (e.g., [1; 2; 3] for a three
%               state model). (Kx1)
%   InitialProbability: Probabality of being in a state in StateSpace
%                       (e.g., [0.4; 0.3; 0.3] for a three state model).
%                       (Kx1)
%   TransitionMatrix: The matrix containing the transition probabilities
%                     from one state to another in a single time step 
%                     (e.g., TransitionMatrix(i, j) = 0.2 means that the
%                     probability of transition from state i to j during a
%                     given time step is 0.2).
%                     Alternatively, TransitionMatrix can have a third
%                     dimension corresponding to the second dimension of
%                     the EmissionMatrix, i.e. the TransitionMatrix is
%                     different for each observation (as might be the case
%                     if the observations were not uniformly sampled in
%                     time).  Note that in this case there will not be a
%                     transition matrix for the first observation.
%                     (numel(StateSpace) x numel(StateSpace), OR
%                      numel(StateSpace) x numel(StateSpace) 
%                          x numel(Observations)-1)
%   EmissionMatrix: The probability density of being in a given state in 
%                   StateSpace given Observation.
%                   (numel(StateSpace) x numel(Observations)).
% 
% OUTPUTS:
%   StateSequence: The estimated state sequence of hidden states given the
%                  input Observations. 

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Ensure inputs are given as column vectors.
if (size(StateSpace, 2) > size(StateSpace, 1))
    StateSpace = StateSpace.';
end
if (size(InitialProbability, 2) > size(InitialProbability, 1))
    InitialProbability = InitialProbability.';
end

% Initialize/pre-allocate variables as needed.
MStates = numel(StateSpace);
NObservations = size(EmissionMatrix, 2);
StateSequence = zeros(NObservations, 1, 'uint8');
StateIndices = StateSequence;
LogPathProbability = zeros(MStates, NObservations);
PathHistory = LogPathProbability;

% Add a third dimension to the TransitionMatrix (if needed) so that code
% doesn't have to be modified with if/else conditions further down.
if (size(TransitionMatrix, 3) == 1)
    TransitionMatrix = repmat(TransitionMatrix, [1, 1, NObservations-1]);
end

% Initialize the log probability table.
LogPathProbability(:, 1) = log(InitialProbability .* EmissionMatrix(:, 1));

% Loop through each observation in Observation and construct the
% log probability table and the path history table.
for jj = 2:NObservations
    for ii = 1:MStates
        [LogPathProbability(ii, jj), PathHistory(ii, jj)] = ...
            max(LogPathProbability(:, jj-1) ...
            + log(TransitionMatrix(:, ii, jj-1)) ...
            + log(EmissionMatrix(ii, jj)));
    end
end

% Construct the state sequence.
[~, StateIndices(NObservations)] = ...
    max(LogPathProbability(:, NObservations));
StateSequence(NObservations) = StateSpace(StateIndices(NObservations));
for ii = NObservations:-1:2
    StateIndices(ii-1) = PathHistory(StateIndices(ii), ii);
    StateSequence(ii-1) = StateSpace(StateIndices(ii-1));
end


end