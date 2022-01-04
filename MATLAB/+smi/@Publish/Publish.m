classdef Publish < handle
    % Publish contains methods useful for batch processing of SR data.
    %   This class contains a collection of analysis/visualization methods
    %   useful for batch-processing of super-resolution data, particularly
    %   for data with multiple labels stored in independent files (e.g.,
    %   sequential super-resolution of two labels on multiple cells).
    %
    % NOTE: This class is designed around .h5 files containing raw data
    %       stored in directories separating distinct cells and labels,
    %       with the directory names following the scheme
    %       Cell*\Label*\Data*.h5
    %
    % REQUIRES: MATLAB 2018a or later (for Publish.genAlignStats())
    %           Image Processing Toolbox
    %           Statistics and Machine Learning Toolbox
    %           Curve Fitting Toolbox
    %
    % CITATION:
    
    % Created by:
    %   David J. Schodt (Lidke Lab 2021), originally based on the script
    %       PublishSEQSR_SAC.m in SR_demo
    
    
    properties
        % Structure of parameters (see smi_core.SingleMoleculeFitting)
        SMF
        
        % Directory containing the Cell*\Label*\Data*.h5 sub-directories.
        CoverslipDir
        
        % Base directory for saving (Default set in performFullAnalysis())
        SaveBaseDir
        
        % Log file for errors (Default set in performFullAnalysis())
        LogFilePath
        
        % Label(s) to be analyzed (Default = [], analyze all labels)
        LabelID = [];
        
        % Cell(s) to be analyzed (Default = [], analyze all cells)
        CellList = [];
        
        % Zoom factor for output SR images (Default = 20)
        SRImageZoom = 20;
        
        % Zoom factor for circle images (Default = 50)
        SRCircleImageZoom = 50;
        
        % Flag to indicate SR results should be generated (Default = true)
        GenerateSR = true;
        
        % Flag to generate various imaging stats (Default = true)
        GenerateImagingStats = true;
        
        % Flag to generate overlay info. between channels (Default = false)
        GenerateOverlayStats = false;
        
        % Flag to perform analysis on bleaching results (Default = false)
        AnalyzeBleaching = false;
        
        % Shift localizations based on brightfield results (Default = true)
        ShiftToReg = true;
        
        % Verbosity of the main analysis workflow. (Default = 1)
        Verbose = 1;
    end
    
    properties (SetAccess = 'protected', Hidden)
        % Instance of SMLM class (for internal use).
        SMLM
        
        % Log of errors encountered during analysis.
        ErrorLog = {};
    end
    
    methods
        function obj = Publish(SMF, CoverslipDir)
            %Publish is the class constructor for the smi.Publish class.
            
            % Set defaults if needed.
            if (~exist('SMF', 'var') || isempty(SMF))
                SMF = smi_core.SingleMoleculeFitting;
            end
            if (~exist('CoverslipDir', 'var') || isempty(CoverslipDir))
                CoverslipDir = '';
            end
            
            % Store the inputs as class properties.
            obj.SMF = SMF;
            obj.CoverslipDir = CoverslipDir;
        end
        
        function set.SMF(obj, SMFInput)
            %set method for the property SMF.
            obj.SMF = smi_core.SingleMoleculeFitting.reloadSMF(SMFInput);
        end
                
        [AlignResultsStruct] = genAlignResults(obj, FilePath, SaveDir);
        genOverlayResults(obj)
        performFullAnalysis(obj)
        processCell(obj, CellName)
        processLabel(obj, CellName, LabelName)
        
    end
    
    methods (Static)
        genSROverlays(ResultsCellDir, SaveDir, AnalysisID)
        [OverlayImage, ColorOrderTag] = overlayNImages(ImageStack);
        genOverlayPlots(ImageShift, RegError, MaxCorr, BPPixelSize, SaveDir)
        [ImagesStruct] = genAlignMovies(AlignRegData, SaveDir);
        [StatsStruct] = genAlignStats(AlignRegStruct, SMD, SaveDir);
        [XCorrStruct] = genAlignXCorr(AlignRegStruct, SaveDir);
        makeOverlayPlots(ImageShift, RegError, MaxCorr, ...
            SRPixelSize, BPPixelSize, SaveDir)
        [PlotAxes, RegError] = plotXYRegError(PlotAxes, SMD);
        [PixelOffsets, SubPixelOffsets, ImageROIs, ImageStats] = ...
            estimateLocalImShifts(Image1, Image2, SubROISize, CorrParams);
        [SubPixelOffsets, SMDROIs, SMDStats] = ...
            estimateLocalCoordShifts(SMD1, SMD2, SubROISize);
        [RegCorrection] = computeRegCorrection(SMF);
        [SMD, BestRegInd] = shiftToBestReg(SMD, RefImage, FocusImages)
    end
    
    
end