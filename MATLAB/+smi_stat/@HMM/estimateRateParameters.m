function [RateParameters, RateParametersSE, LogLikelihood] = ...
    estimateRateParameters(EmissionPDFCell, DeltaT, ...
    RateParametersGuess, SearchOptions, Verbose)
%estimateRateParameters estimates rate parameters from dimer candidates.
% This method will perform a Hidden Markov Model (HMM) analysis on all
% trajectory pairs contained in TRArray, using the model given by the 
% emission probability densities EmissionPDFCell.
% 
% INPUTS:
%   EmissionPDFCell: A cell array containing the emission PDF's for each
%                    state for every candidate in TRArray.  For example,
%                    EmissionPDFCell{10} will be equal to
%                    TRArray(1, 10).EmissionProbabilities(...
%                       TRArray(1, 10).DimerCandidateBool)
%   DeltaT: Cell array with each element corresponding to the time since
%           last observation, corresponding to the observations used to 
%           generate EmissionPDFCell (e.g., DeltaT{n} = diff(TRArray(1,
%           n).FrameNum(TRArray(1, n).DimerCandidateBool)).
%           (NCandidates x 1 cell array)
%   MaxSeparation: Maximum separation that may exist between trajectories
%                  in TRArray.
%   RateParametersGuess: An initial guess of the rate parameters. 
%                        (default = 0.1 * ones(NRateParameters, 1))
%                        (NRateParameters x 1 array, where 
%                         NRateParameters = NStates * (NStates-1))
%   SearchOptions: MATLAB optimoptions structure sent directly to
%                  fmincon(). (Default = optimoptions('fmincon'))
%   Verbose: Verbosity level requested. (Default = 1)
%
% OUTPUTS:
%   RateParameters: An array containing the rate parameters estimated for
%                   the HMM (see the nested function 
%                   computeLogLikelihoodWrapper below for parameter
%                   descriptions).
%                   (1 / frames)
%   RateParametersSE: The estimated standard errors in the reported 
%                     RateParameters. (1 / frames)
%   LogLikelihood: The log-likelihood computed at the returned value of
%                  RateParameters. 
%
% REQUIRES: Optimization Toolbox to use fmincon()

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Define misc. parameters/set defaults.
if (~exist('Verbose', 'var') || isempty(Verbose))
    Verbose = 1;
end
if (~exist('SearchOptions', 'var') || isempty(SearchOptions))
    SearchOptions = optimoptions('fmincon');
    SearchOptions.Display = ...
        smi_helpers.arrayMUX({'none', 'final', 'iter'}, Verbose);
end
NStates = size(EmissionPDFCell{1}, 2);
NRateParams = NStates * (NStates-1);

% Set defaults if needed.
if ~exist('RateParametersGuess', 'var') || isempty(RateParametersGuess)
    RateParametersGuess = 0.01 * ones(NRateParams, 1);
end

% Estimate the rate parameters by minimizing the negative of the
% log-likelihood.
NegLogLikelihood = @(X) -computeLogLikelihoodWrapper(X, DeltaT, ...
    EmissionPDFCell, NStates);
RateParameters = fmincon(NegLogLikelihood, RateParametersGuess, ...
    [], [], [], [], zeros(NRateParams, 1), inf(NRateParams, 1), ...
    [], SearchOptions);
LogLikelihood = -NegLogLikelihood(RateParameters);

% Estimate the errors in the found model parameters.
[HessianMatrix] = smi_stat.computeHessian(NegLogLikelihood, RateParameters);
RateParametersSE = sqrt(diag(inv(HessianMatrix)));


    function [LogLikelihood] = computeLogLikelihoodWrapper(...
            SearchParameters, DeltaT, EmissionPDFCell, NStates)
        % This is a wrapper for the class method computeLogLikelihood.  In
        % computeLogLikelihood, the rate parameters are given as a matrix,
        % but to estimate the rate parameters with fminsearch we need to
        % have a function with a vector input.  This function also acts to
        % loop over all trajectory pairs to compute the overall likelihood
        % of a given model parametrized by RateParameters.
        % 
        % INPUTS: 
        %   SearchParameters: A vector containing the transition rate
        %                     parameters between the HMM states as well as
        %                     the domain size parameter.
        %                     ModelSpecifier = 'DF': 
        %                       SearchParameters(1) = rate parameter from
        %                           state 2 to state 1 (Free to Dimer)
        %                       SearchParameters(2) = rate parameter from
        %                           state 1 to state 2 (Dimer to Free)
        %                     ModelSpecifier = 'DDF':
        %                       SearchParameters(1) = rate parameter from
        %                           state 2 to state 1 (Domain to Dimer)
        %                       SearchParameters(2) = rate parameter from
        %                           state 3 to state 1 (Free to Dimer)
        %                       SearchParameters(3) = rate parameter from
        %                           state 1 to state 2 (Dimer to Domain)
        %                       SearchParameters(4) = rate parameter from
        %                           state 3 to state 2 (Free to Domain)
        %                       SearchParameters(5) = rate parameter from
        %                           state 1 to state 3 (Dimer to Free)
        %                       SearchParameters(6) = rate parameter from
        %                           state 2 to state 3 (Domain to Free)
        %                   NOTE: The ordering of rate parameters is based
        %                         on MATLABs column based behavior, so to 
        %                         convert from matrix to array we scan down 
        %                         each column and read off the elements.
        %   EmissionPDFCell: A cell array containing the emission PDF's for
        %                    each state for all trajectory pairs in TRArray
        %   NStates: Number of states in the model.
        %   
        % OUTPUTS:
        %   LogLikelihood: The log-likelihood of the sequence as computed
        %                  by smi_stat.HMM.computeLogLikelihood()
                
        % Convert the SearchParameters array to a rate parameter matrix for
        % use in smi_stat.HMM.computeLogLikelihood(). 
        RateParameterMatrix = zeros(NStates);
        DiagonalBool = logical(eye(NStates));
        RateParameterMatrix(~DiagonalBool) = SearchParameters;
                
        % Loop through all dimer candidate pairs in TRArray to compute the
        % log-likelihood.
        LogLikelihood = 0;
        for nn = 1:numel(DeltaT)
            % Compute the transition matrices for all observations of the 
            % current trajectory pair.
            TransitionMatrixSeries = smi_stat.HMM.generateTransitionMatrix(...
                RateParameterMatrix, DeltaT{nn});
                        
            % Compute the log-likelihood.
            LogLikelihood = LogLikelihood ...
                + smi_stat.HMM.computeLogLikelihood(...
                TransitionMatrixSeries, EmissionPDFCell{nn});
        end
    end


end