classdef SMLM < handle
    % Single Molecule Localization Microscopy Analysis
    %
    % This is a high-level class that provides complete analysis of SMLM data.
    %
    %
    
    % =========================================================================
    properties
        SMF
        Preset      % {'TIRF', 'Sequential'} good idea?
        Data        % Current dataset or used for manual setting of data
        DataType    % {'File', 'UserDefined'} ?
        FileName    % String or Cell array of strings
        DataDir
        ResultsDir  % (Default = 'DataDir/../Results/FileName/) same as Seq
    end
    % =========================================================================
    
    % =========================================================================
    properties (Access=protected)
        DC  % DriftCorrection class object used internally
        SMD % SMD structure with final analysis results
        NDataSets %
    end % properties (Access=protected)
    % =========================================================================
    
    % =========================================================================
    methods
        function obj=SMLM(SMF,Filename)
            % SMLM
        end
        
        function fullAnalysis(obj)
            % fullAnalysis Analyze all data and save results
            
            %analyzeAll
            
            %saveResults
            
            
            
            %save
            
        end
        
        function analyzeAll(obj)
            % analyzeAll loops over dataset and creates SMD
            
            % DriftCorrection class object is also used in analyzeDataset
            obj.DC = smi_core.DriftCorrection(obj.SMF);
            obj.SMD=[];
            for nn=1:obj.NDatasets
                SMDnn=obj.analyzeDataset(nn);
                obj.SMD=smi_core.SingleMoleculeData.catSMD(obj.SMD,SMDnn);
            end
            
            % Inter-dataset drift correction.
            fprintf('Drift correcting (inter-datastet) ...\n');
            obj.SMD = obj.DC.driftCorrectKNNInter(obj.SMD);
        end
        
        
        function SMD=analyzeDataset(obj,DataSetIndex)
        % analyzeDataset Load and analyze one dataset    
            
            obj.SMF=smi_core.SingleMoleculeFitting.createSMF();
            fprintf('Loading dataset %d ...\n', DataSetIndex);
            [Dataset, obj.SMF]=obj.loadDataset(obj.SMF,DataSetIndex);
            
            % Generate localizations from the current Dataset.
            LD = smi_core.LocalizeData(Dataset, obj.SMF);
            fprintf('Generating localizations ...\n');
            [SMD] = LD.genLocalizations();
            
            % Perform frame-connection on localizations in SMD.
            if obj.SMF.FrameConnection.On 
                FC = smi_core.FrameConnection(SMD, obj.SMF);
                fprintf('Frame connecting ...\n');
                [SMD] = FC.performFrameConnection();
            end
            
            % Intra-dataset drift correction.
            if obj.SMF.DriftCorrection.On 
                fprintf('Drift correcting (intra-datastet) ...\n');
                SMD = obj.DC.driftCorrectKNNIntra(SMD, DataSetIndex);
            end
        end
        
        
        function [Dataset, SMF]=loadDataset(obj,SMF,DataSetIndex)
        % loadDataset loads a dataset and converts to photons
        % set obj.Data   
        [~, Dataset, SMF] = smi_core.LoadData(SMF,DataSetIndex);
        end
        
        function saveResults(obj)
        % saveResults Save all results and plots in subfolder
        % gaussblobs, drift image, fits/frame, NumConnected hist,
        % Driftcorrection plots, precision hist, intensity hist,
        % mat file with SMD and SMF files. 
        
        
        end
        
    end % methods
    % =========================================================================


    % =========================================================================
    methods(Static)
        Success = unitTest();
    end % methods(Static)
    % =========================================================================
    
end
