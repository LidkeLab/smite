function [RateParameters, RateParametersSE, LogLikelihood] = ...
    performFullAnalysis(obj)
%performFullAnalysis performs all analyses on the data in obj.TRArray.
% This method performs all analyses possible on the data.  This is intended
% to be a main "run" method to this class, meaning that the user can
% create an instance of this class, set all class properties as needed, and
% then run this method to perform a complete analysis of the data.
% NOTE: Outputs from this method are added for convenience:
%       they should all also be saved/updated in the class instance.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Add a message to the Command Window to indicate the status.
if (obj.Verbose > 0)
    fprintf('HMM.performFullAnalysis(): beginning analysis...\n')
end

% Make local copies of various arrays and define misc. parameters.
TRArray = obj.TRArray;
NStates = numel(obj.PDFHandles);
if ~NStates
    error('You must set smi_stat.HMM.PDFHandles before proceeding')
end
StateSpace = (1:NStates).';

% Isolate the dimer candidates that we wish to send through the HMM
% analysis and add some new fields that we'll need.
TRArrayTrunc = obj.isolateCandidateTRArray(TRArray);
NCandidates = size(TRArrayTrunc, 1);

% Ensure that the arrays obj.DiffusionCoefficient and obj.RegistrationError
% make sense based on the size of TRArray.
NRegError = numel(obj.RegistrationError);
if (NRegError ~= NCandidates)
    % In this case, we'll take the first entry of obj.RegistrationError and
    % use it for all candidates.
    obj.RegistrationError = ...
        obj.RegistrationError(1) * ones(NCandidates, 1);
end
DSize = size(obj.DiffusionCoefficient);
if ((DSize(1)==NCandidates) && (DSize(2)==2))
    % One diffusion coefficient was provided per trajectory, so no action
    % is needed.
elseif (prod(DSize) == 2)
    % If two diffusion coefficients were provided, we'll use those values
    % for each trajectory (depending on the channel).
    obj.DiffusionCoefficient = ...
        [obj.DiffusionCoefficient(1) * ones(NCandidates, 1), ...
        obj.DiffusionCoefficient(2) * ones(NCandidates, 1)];
else
    % In this case, we'll just take the first element of
    % DiffusionCoefficient and use it for all trajectories.
    obj.DiffusionCoefficient = ...
        obj.DiffusionCoefficient(1) * ones(NCandidates, 2);
end

% Compute the emission pdf's for the HMM states for all of the trajectory
% pairs.
if (obj.Verbose > 1)
    fprintf('HMM.performFullAnalysis(): preparing emission matrix...\n')
end
EmissionPDFCell = cell(NCandidates, 1);
EmissionPDFInputs = cell(8, 1);
EmissionPDFInputs{5} = obj.DimerSeparation;
EmissionPDFInputs{7} = obj.MaxSeparation;
EmissionPDFInputs{8} = obj.DomainSeparation;
DeltaT = EmissionPDFCell;
for ii = 1:NCandidates
    % Update the EmissionPDFInputs for this candidate.
    EmissionPDFInputs{1} = TRArrayTrunc(ii, 1).Separation;
    EmissionPDFInputs{2} = [TRArrayTrunc(ii, 1).AverageSE, ...
        TRArrayTrunc(ii, 2).AverageSE];
    EmissionPDFInputs{3} = double(diff(TRArrayTrunc(ii, 1).FrameNum));
    EmissionPDFInputs{4} = obj.RegistrationError(ii);
    EmissionPDFInputs{6} = obj.DiffusionCoefficient(ii, :);
    DeltaT{ii} = EmissionPDFInputs{3};
    
    % Compute the emission pdf's and store them in the TRArray.
    EmissionPDFCell{ii} = ...
        obj.generateEmissionMatrix(obj.PDFHandles, EmissionPDFInputs);
    TRArrayTrunc(ii, 1).EmissionProbabilities = EmissionPDFCell{ii};
    TRArrayTrunc(ii, 2).EmissionProbabilities = EmissionPDFCell{ii};
end

% Estimate the rate parameters from the data and convert to a matrix for
% later use.
if (obj.Verbose > 1)
    fprintf('HMM.performFullAnalysis: estimating rate parameters...\n')
end
[RateParameters, RateParametersSE, LogLikelihood] = ...
    obj.estimateRateParameters(EmissionPDFCell, DeltaT, ...
    obj.RateParametersGuess);
TransitionRates = zeros(NStates);
DiagonalBool = logical(eye(NStates));
TransitionRates(~DiagonalBool) = RateParameters;

% Estimate the state sequence for each pair of trajectories, storing the
% result in TRArray.
if (obj.Verbose > 1)
    fprintf('HMM.performFullAnalysis: computing state sequence...\n')
end
for ii = 1:NCandidates
    % Generate the transition matrices for these trajectories.
    TransitionMatrixSeries = obj.generateTransitionMatrix(...
        TransitionRates, diff(TRArrayTrunc(ii, 1).FrameNum));
    
    % Estimate the state sequence that resulted in these observations.
    EmissionProbabilitySeries = EmissionPDFCell{ii}.';
    [StateSequence] = obj.estimateStateSequence(StateSpace, ...
        EmissionProbabilitySeries(:, 1), TransitionMatrixSeries, ...
        EmissionProbabilitySeries);
    
    % Find the apparent durations of each dimer event (as found by the
    % Viterbi algorithm in obj.estimateStateSequence() above).
    DimerDurations = obj.computeDimerDurations(...
        StateSequence, TRArrayTrunc(ii, 1).FrameNum);
    
    % Store interesting arrays in the TRArrayTrunc.
    TRArrayTrunc(ii, 1).StateSequence = StateSequence;
    TRArrayTrunc(ii, 1).DimerDurations = DimerDurations;
    TRArrayTrunc(ii, 2).StateSequence = StateSequence;
    TRArrayTrunc(ii, 2).DimerDurations = DimerDurations;
    
    % Store interesting arrays in the "full sized" TRArray (useful for
    % some plots/movies we might wish to make).
    DimerCandidateBool1 = TRArray(ii, 1).DimerCandidateBool;
    DimerCandidateBool2 = TRArray(ii, 2).DimerCandidateBool;
    Separation1Padded = NaN(numel(DimerCandidateBool1), 1);
    Separation2Padded = NaN(numel(DimerCandidateBool2), 1);
    Separation1Padded(DimerCandidateBool1) = ...
        TRArrayTrunc(ii, 1).Separation;
    Separation2Padded(DimerCandidateBool2) = ...
        TRArrayTrunc(ii, 2).Separation;
    StateSequence1Padded = NaN(numel(DimerCandidateBool1), 1);
    StateSequence2Padded = NaN(numel(DimerCandidateBool2), 1);
    StateSequence1Padded(DimerCandidateBool1) = StateSequence;
    StateSequence2Padded(DimerCandidateBool2) = StateSequence;
    EmissionProb1Padded = NaN(numel(DimerCandidateBool1), NStates);
    EmissionProb2Padded = NaN(numel(DimerCandidateBool2), NStates);
    EmissionProb1Padded(DimerCandidateBool1, :) = ...
        TRArrayTrunc(ii, 1).EmissionProbabilities;
    EmissionProb2Padded(DimerCandidateBool2, :) = ...
        TRArrayTrunc(ii, 2).EmissionProbabilities;
    TRArray(ii, 1).Separation = Separation1Padded;
    TRArray(ii, 1).StateSequence = StateSequence1Padded;
    TRArray(ii, 1).DimerDurations = DimerDurations;
    TRArray(ii, 1).EmissionProbabilities = EmissionProb1Padded;
    TRArray(ii, 2).Separation = Separation2Padded;
    TRArray(ii, 2).StateSequence = StateSequence2Padded;
    TRArray(ii, 2).DimerDurations = DimerDurations;
    TRArray(ii, 2).EmissionProbabilities = EmissionProb2Padded;
end

% Determine if the DimerDurations and the found dimer off-rate are in
% agreement.
DimerDurations = cell2mat({TRArrayTrunc(:, 1).DimerDurations}.');
ViterbiOffRate = -log(1 - 1/mean(DimerDurations));
if (NStates == 2)
    HMMOffRate = RateParameters(2);
elseif (NStates == 3)
    HMMOffRate = RateParameters(3) + RateParameters(5);
else
    error('Update performFullAnalysis() for new state model!!!');
end
RelativeDiff = (abs(ViterbiOffRate-HMMOffRate) ...
    / abs(ViterbiOffRate+HMMOffRate));
if ((RelativeDiff>obj.DiscrepancyTol) && (obj.Verbose>0))
    warning(['Off-rate estimates differ by %0.2g%%, ', ...
        '%0.2g%% tolerance not met'], ...
        100*RelativeDiff, 100*obj.DiscrepancyTol);
end

% Save misc. local arrays back into the class instance.
obj.RateParameters = RateParameters;
obj.RateParametersSE = RateParametersSE;
obj.TRArray = TRArray;
obj.TRArrayTrunc = TRArrayTrunc;

% Save the analysis results.
if obj.SaveResults
    if (obj.Verbose > 1)
        fprintf('HMM.performFullAnalysis: saving results...\n')
    end
    obj.saveResults();
end

% Add an update message to the Command Window.
if (obj.Verbose > 0)
    fprintf('HMM.performFullAnalysis(): analysis complete!\n')
end


end