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
        
        CoverslipDir
        SaveBaseDir
        LabelID = []; % label(s) to analyze, if empty analyze all labels
        GenerateSR = 1; % flag to generate SR images in the analysis
        GenerateImagingStats = 1; % flag to generate misc. imaging info.
        GenerateOverlayStats = 0; % flag to generate misc. overlay info.
        AnalyzeBleaching = 0; % flag to analyze bleaching acquisitions

        
        % Verbosity of the main analysis workflow. (Default = 1)
        Verbose = 1;
    end
    
    properties (SetAccess = protected)
        PublishedResultsStruct = struct(); % struct. w/ published results
        CellLabelStruct = struct(); % struct. w/ concat. published results
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
        
        performFullAnalysis(obj)
        processCell(obj)
        processLabel(obj)
        saveResults(obj)
        
    end
    
    methods(Static)
        
    end
    
    
end