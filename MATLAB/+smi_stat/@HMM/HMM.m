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
    % & Lidke, D. S. (2021). EGFR transactivates RON to drive oncogenic
    % crosstalk. eLife
    % https://doi.org/10.7554/eLife.63678
    %
    % REQUIRES:
    %   Optimization Toolbox
    %   Image Processing Toolbox
    %   Statistics and Machine Learning Toolbox
    
    % Created by:
    %   David J. Schodt (Lidke Lab, 2021)
    
    
    properties
        % Separation between fluorophores on a dimer (pixels)
        DimerSeparation(1, 1) {mustBeFloat(DimerSeparation)} = 0.5;
        
        % Typical domain size for the free, dimer, domain model (pixels)
        DomainSeparation(1, 1) {mustBeFloat(DomainSeparation)} = 2;
        
        % Max. separation for dimer candidates (pre-processing) (pixels)
        MaxSeparation(1, 1) {mustBeFloat(MaxSeparation)} = 5;

        % Registration error inflation (pixels)
        RegErrorInflation(1, 1) {mustBeFloat(RegErrorInflation)} = 0.0;
        
        % Identifier for one of the pre-built models. (Default = 'DF')
        % OPTIONS:
        %   'DF': dimer or free
        %   'DDF': dimer, domain, or free
        ModelSpecifier {mustBeText(ModelSpecifier)} = 'DF';
        
        % Handles to the state PDFs used in the HMM (cell array)
        PDFHandles(:, 1) cell
        
        % Initial guess of rate parameters (NRatesx1 float)
        RateParametersGuess(:, 1) {mustBeFloat(RateParametersGuess)} = ...
            0.01 * [1; 1];
        
        % Array of TR structures corresponding to dimer candidates. (Nx2)
        % This is organized as a NCandidatex2 structure, with each column
        % being a dimer candidate.
        TRArray(:, 2) struct
        
        % Data channel names added to certain outputs. (cell array of char)
        ChannelNames cell = {'Channel 1'; 'Channel 2'};
        
        % Model state names added to certain outputs. (cell array of char)
        StateNames cell = {'Dimer'; 'Domain'; 'Free'};
        
        % Label for save directory to indicate a condition. (char/string)
        ConditionLabel {mustBeText(ConditionLabel)} = '';
        
        % Indicate results should be saved. (Default = true)
        % NOTE: This is used when running obj.performFullAnalysis().
        SaveResults logical = true;
        
        % Indicate outputs should be in physical units. (Default = false)
        UnitFlag logical = false;
        
        % Indicate movies should be generated and saved. (Default = false)
        % NOTE: This is used when running obj.performFullAnalysis().
        GenerateMovies logical = false;
        
        % Set of movie parameters (see GenerateMovies)
        MovieParams struct = smi_vis.GenerateMovies.prepDefaults();
        
        % Set of plot parameters (see DisplayParams in plotDimerPairInfo())
        PlotParams struct = struct();
        
        % Indicate plots should be generated and saved. (Default = true)
        % GeneratePlots(1) indicates the histogram of dimer durations
        %   should be generated.
        % GeneratePlots(2) indicates the summaray plots should be
        %   generated.
        % NOTE: This is used when running obj.performFullAnalysis() inside
        %       of obj.saveResults().
        GeneratePlots = [true; true];

        % Structure of parameters (see smi_core.SingleMoleculeFitting)
        % NOTE: As of this writing, this class uses SMF.Data.FrameRate and
        %       SMF.Data.PixelSize (when UnitFlag = true).  If
        %       DiffusionCoefficient is empty, obj.performFullAnalysis()
        %       will use the value in SMF.Tracking.D.
        SMF
        
        % Top level directory for saving results.
        % NOTE: If left empty, obj.saveResults() will try to use
        %       obj.SMF.Data.ResultsDir.  If that is also empty, a default
        %       will be set to pwd().
        SaveDir {mustBeText(SaveDir)} = '';
        
        % Verbosity level of obj.performFullAnalysis(). (Default = 1)
        Verbose {mustBeInteger(Verbose)} = 1;
    end
    
    properties (SetAccess = 'protected')
        % Rate parameters found by HMM analysis. (float array)
        RateParameters(:, 1) {mustBeFloat(RateParameters)}
        
        % Standard error estimates of rate parameters. (float array)
        RateParametersSE(:, 1) {mustBeFloat(RateParametersSE)}
        
        % Pre-processed obj.TRArray as seen by the HMM analysis. (Nx2)
        TRArrayTrunc(:, 2) struct
    end
    
    properties (Hidden)
        % Tolerance used in a check in obj.performFullAnalysis().
        DiscrepancyTol(1, 1) {mustBeFloat(DiscrepancyTol)} = 0.1;
        
        % Filenames corresponding to dimer candidates in TRArray.
        FileNames(:, 2) cell
        
        % Registration files corresponding to dimer candidates in TRArray.
        RegFileNames(:, 1) cell
    end
    
    methods
        function obj = HMM(TRArray, SMF)
            % Class constructor which allows option inputs.
            if exist('TRArray', 'var')
                obj.TRArray = TRArray;
            end
            if exist('SMF', 'var')
                obj.SMF = SMF;
            end

            % Revise some default movie parameters.
            obj.MovieParams.AutoCrop = true;
            obj.MovieParams.IndicateDimer = true;
            obj.MovieParams.IndicateDimerCandidate = true;
        end
        
        function set.ConditionLabel(obj, InputValue)
            % This set method ensures that the user defined ConditionLabel
            % doesn't contain any spaces (spaces can cause issues in some
            % parts of the code, but are a common choice by users, so I've
            % found it easiest to just replace them quietly).
            obj.ConditionLabel = replace(InputValue, ' ', '_');
        end
        
        [RateParameters, RateParametersSE, LogLikelihood] = ...
            performFullAnalysis(obj);
        saveResults(obj);
        
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
        [TRArray, FileList] = ...
            findDimerCandidatesFromFiles(FileDir, FilePatterns, ...
            MaxDimerSeparation, MaxSeparation, MinValidPoints, ...
            MinPhotons, BorderPadding, Verbose);
        [TRArray] = computePairSeparation(TRArray);
        [TRArrayTrunc] = isolateCandidateTRArray(TRArray);
        [StateSequence] = computeViterbiPath(...
            StateSpace, InitialProbability, ...
            TransitionMatrix, EmissionMatrix);
        [DimerDurations] = computeDimerDurations(...
            StateSequence, FrameNum, DimerState);
        [PlotAxes] = plotOffRateBarGraph(OffRates, OffRatesSE, ...
            XBarLocations, ConditionColorMap, ConditionNames, UnitFlag, ...
            PlotAxes);
        [PlotAxes] = plotDimerDurationCDF(DimerDurations, PlotAxes);
        [MovieParams] = createDimerMovie(MovieAxes, ...
            TRArray, RawDataChannel1, RawDataChannel2, ...
            FilePath, MovieParams, SMF, VideoObject)
        [MovieParams] = createAllMovies(TRArray, ...
            MovieParams, SaveDir, RawDataBaseDir);
        [FigureHandle, DisplayParams] = createSummaryPlot(FigureHandle, ...
            TRArray, SMF, DisplayParams, UnitFlag);
    end
    
    methods (Static, Hidden)
        % These methods should still be accesible to the user, but it's
        % probably best that we don't overwhelm the user with too many
        % visible methods that they aren't likely to use!
        
        [PlotFigure, DisplayParams] = plotSepDistribs(PlotFigure, ...
            TRArray, PDFHandles, PDFInputs, DisplayParams);
        [PlotAxes, DisplayParams] = plotDimerPairInfo(PlotAxes, ...
            PlotType, TRArray, SMF, DisplayParams, UnitFlag);
        %         plotViterbiPath
        %         plotEmissionProbabilities
        %         plotDimerTraj2D
        %         plotDimerTraj3D
        %         plotRegistration
        %         plotXYSeparation
    end
    
    
end