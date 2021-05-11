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
    % REQUIRES: MIC_H5 class from the matlab-instrument-control repository
    %           MATLAB 2018a or later (for Publish.genAlignStats())
    %           MATLAB Image Processing Toolbox 2014a or later
    %           DipImage (to use findshift())
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
        
        % Label(s) to be analyzed (Default = [], analyze all labels)
        LabelID = [];
        
        % Zoom factor for output SR images (Default = 20)
        SRImageZoom = 20;
        
        % Flag to indicate SR results should be generated (Default = true)
        GenerateSR = true;
        
        % Flag to generate various imaging stats (Default = true)
        GenerateImagingStats = true;
        
        % Flag to generate overlay info. between channels (Default = false)
        GenerateOverlayStats = false;
        
        % Flag to perform analysis on bleaching results (Default = false)
        AnalyzeBleaching = false;
        
        % Verbosity of the main analysis workflow. (Default = 1)
        Verbose = 1;
        
        % Structure containing several analysis results.
        PublishedResultsStruct = struct();
        
        % Structure containing several concatenated analysis results.
        CellLabelStruct = struct();
    end
    
    properties (SetAccess = 'protected', Hidden)
        % Instance of SMLM class (for internal use).
        SMLM
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
        
        function saveResults(obj)
            % This method saves the PublishedResultsStruct into the Results
            % directory.
            PublishedResultsStruct = obj.PublishedResultsStruct;
            CellLabelStruct = obj.CellLabelStruct;
            save(fullfile(obj.SaveBaseDir, 'ResultsStruct.mat'), ...
                'PublishedResultsStruct');
            save(fullfile(...
                obj.SaveBaseDir, 'ResultsStructConcatenated.mat'), ...
                'CellLabelStruct');
        end
        
        [AlignResultsStruct] = genAlignResults(obj, FilePath, SaveDir);
        performFullAnalysis(obj)
        processCell(obj, CellName)
        processLabel(obj, CellName, LabelName)
        
    end
    
    methods(Static)
        genSROverlays(ResultsCellDir, SaveDir)
        [OverlayImage, ColorOrderTag] = overlayNImages(ImageStack);
        genOverlayPlots(ImageShift, RegError, MaxCorr, SRPixelSize, ...
            BPPixelSize, SaveDir)
        [ImagesStruct] = genAlignMovies(AlignRegData, SaveDir);
        [StatsStruct] = genAlignStats(AlignRegStruct, SMDR, SaveDir);
        [XCorrStruct] = genAlignXCorr(AlignRegStruct, SaveDir);
        makeOverlayPlots(ImageShift, RegError, MaxCorr, ...
            SRPixelSize, BPPixelSize, SaveDir)
        [CellLabelStruct] = concatenateResults(PublishedResultsStruct);
        genConcatenatedFigures(CellLabelStruct, SaveDir);
        [PlotAxes, RegError] = plotXYRegError(PlotAxes, SMD);
        [PixelOffsets, SubPixelOffsets, BlockROIs] = ...
            estimateLocalShifts(Image1, Image2, BlockSize, UseGPU);
    end
    
    
end