classdef HMM < handle
    %HMM contains methods for hidden Markov model analysis.
    % This class contains methods related to/useful for hidden Markov model
    % (HMM) of single-particle tracking (SPT).  In particular, this class
    % is designed to perform HMM and related analyses on the output
    % Tracking Results (TR) structures produced by smi.SPT.
    %
    % NOTE: This class ALWAYS uses camera units (pixels, frames)
    %       throughout the analysis.  When UnitFlag = 1, units are not
    %       converted to physical units (micrometers, seconds) until the
    %       very end when outputs/plots are being produced.
    %
    % See reference 
    % Nitta, C. F., Green, E. W., Jhamba, E. D., Keth, J. M.,
    % Ortiz-Caraveo, I., Grattan, R. M., Schodt, D. J., Gibson, A. C., 
    % Rajput, A., Lidke, K. A., Steinkamp, M. P., Wilson, B. S., 
    % & Lidke, D. S. (2020). EGFR transactivates RON to drive oncogenic 
    % crosstalk. BioRxiv, 2020.08.11.246785. 
    % https://doi.org/10.1101/2020.08.11.246785
    %
    % REQUIRES: Optimization Toolbox
    
    % Created by:
    %   David J. Schodt (Lidke Lab, 2021)
    
    
    properties
        % Diffusion constant(s) for the trajectories (pixel^2 / frame)
        DiffusionConstant {mustBeFloat(DiffusionConstant)}
        
        % Separation between fluorophores on a dimer (pixels)
        DimerSeparation(1, 1) {mustBeFloat(DimerSeparation)} = 0.5;
        
        % Typical domain size for the free, dimer, domain model (pixels)
        DomainSeparation(1, 1) {mustBeFloat(DomainSeparation)} = 2;
        
        % Max. separation for dimer candidates (pre-processing) (pixels)
        MaxSeparation = 5;
        
        % Handles to the state PDFs used in the HMM (cell array)
        PDFHandles cell
        
        % Initial guess of rate parameters (NRatesx1 float)
        RateParametersGuess {mustBeFloat(RateParametersGuess)}
        
        % Array of TR structures corresponding to dimer candidates (2xN)
        % This is organized as a 2xNCandidate structure, with each column
        % being a dimer candidate.
        TRArray struct
        
        % Data channel names added to certain outputs. (cell array of char)
        ChannelNames cell = {'Channel 1', 'Channel 2'};
        
        % Model state names added to certain outputs. (cell array of char)
        StateNames cell = {'Dimer', 'Domain', 'Free'};
        
        % Label for save directory to indicate a condition. (char/string)
        ConditionLabel {mustBeText(ConditionLabel)}
        
        % Indicate results should be saved. (Default = true)
        % NOTE: This is used when running obj.performFullAnalysis().
        SaveResults logical = true;
        
        % Indicate outputs should be in physical units. (Default = false)
        UnitFlag logical = false;
        
        % Indicate movies should be generated and saved. (Default = false)
        % NOTE: This is used when running obj.performFullAnalysis().
        GenerateMovies logical = false;
        
        % Indicate plots should be generated and saved. (Default = true)
        % NOTE: This is used when running obj.performFullAnalysis().
        GeneratePlots logical = true;
        
        % Top level directory for saving results.
        SaveDir {mustBeText(SaveDir)}
    end
    
    properties (SetAccess = 'protected')
        % Rate parameters found by HMM analysis. (float array)
        RateParameters {mustBeFloat(RateParameters)}
        
        % Standard error estimates of rate parameters. (float array)
        RateParametersSE {mustBeFloat(RateParametersSE)}
        
        % Pre-processed obj.TRArray as seen by the HMM analysis. (2xN)
        TRArrayTrunc struct
    end
    
    methods
        function obj = HMM()
            
        end
        
        [RateParameters, RateParametersSE, LogLikelihood] = ...
            performFullAnalysis(obj);
        saveResults(obj);
        
        function set.ConditionLabel(obj, InputValue)
            % This set method ensures that the user defined ConditionLabel
            % doesn't contain any spaces (spaces can cause issues in some
            % parts of the code, but are a common choice by users, so I've
            % found it easiest to just replace them quietly).
            obj.ConditionLabel = replace(InputValue, ' ', '_');
        end
        
    end
    
    methods (Static)
        [TransitionMatrix] = ...
            generateTransitionMatrix(TransitionRates, DeltaT);
        [EmissionMatrix] = ...
            generateEmissionMatrix(PDFHandles, DataCell);
        [PDFHandles] = generateEmissionPDFs(ModelSpecifier);
        [LogLikelihood] = computeLogLikelihood(TransitionMatrixSeries, ...
            EmissionProbabilitySeries);
        [RateParameters, RateParametersSE, LogLikelihood] = ...
            estimateRateParameters(EmissionPDFCell, DeltaT, ...
            RateParametersGuess, SearchOptions, Verbose);
        [TRArray] = findDimerCandidates(TR1, TR2, ...
            MaxDimerSeparation, MaxSeparation, MinValidPoints, ...
            MinPhotons, BorderPadding);
        [TRArray] = computePairSeparation(TRArray);
        [TRArrayTrunc] = isolateCandidateTRArray(TRArray);
        [StateSequence] = computeViterbiPath(...
            StateSpace, InitialProbability, ...
            TransitionMatrix, EmissionMatrix);
        [DimerDurations] = computeDimerDurations(...
            StateSequence, FrameNum, DimerState);
        [PlotAxes, DisplayParams] = plotDimerPairInfo(TRArray, ...
            DisplayParams, PlotType, PlotAxes);
        [PlotAxes] = plotOffRateBarGraph(OffRates, OffRatesSE, ...
            XBarLocations, ConditionColorMap, ConditionNames, UnitFlag, ...
            PlotAxes);
        [PlotAxes] = plotDimerDurationCDF(DimerDurations, PlotAxes);
        [DisplayParams] = createDimerMovie(TRArray, ...
            RawDataChannel1, RawDataChannel2, ...
            FilePath, DisplayParams, PlotAxes);
        [DisplayParams] = createAllMovies(TRArray, ...
            SaveDir, RawDataBaseDir, DisplayParams);
        [FigureHandle] = createSummaryPlot(TRArray, DisplayParams, ...
            FigureHandle);
    end
    
end