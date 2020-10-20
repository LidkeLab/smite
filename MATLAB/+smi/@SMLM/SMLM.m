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
    end % properties (Access=protected)
    % =========================================================================

    % =========================================================================
    methods
        function obj=SMLM(SMF,Filename)
            % SMLM
            obj.SMF = SMF;
        end
        
        function fullAnalysis(obj)
            % fullAnalysis Analyze all data and save results
            
            obj.analyzeAll();
            obj.saveResults();
            
            %save
            
        end

        % ---------------------------------------------------------------------

        function analyzeAll(obj)
            % analyzeAll loops over dataset and creates SMD.
            
            % Define the list of datasets to be processed.
            %obj.SMF.Data.DatasetList = ...
            %    smi_core.LoadData.setDatasetList(obj.SMF).
            %datasetList = obj.SMF.Data.DatasetList;
            datasetList = [1, 2];

            % DriftCorrection class object is also used in analyzeDataset
            obj.DC = smi_core.DriftCorrection(obj.SMF);
            obj.SMD=[];
            for nn=1:numel(datasetList)
                SMDnn = obj.analyzeDataset(datasetList(nn), nn);
                obj.SMD=smi_core.SingleMoleculeData.catSMD(obj.SMD,SMDnn);
            end
            
            % Inter-dataset drift correction.
            if numel(datasetList) > 1
               fprintf('Drift correcting (inter-datastet) ...\n');
               obj.SMD = obj.DC.driftCorrectKNNInter(obj.SMD);
            end
        end
        
        % ---------------------------------------------------------------------
        
        function SMD=analyzeDataset(obj,DatasetIndex,DatasetCount)
        % analyzeDataset Load and analyze one dataset    

            if ~exist('DatasetCount', 'var')
                DatasetCount = 1;
            end
            
            fprintf('Loading dataset %d ...\n', DatasetIndex);
            [Dataset, obj.SMF]=obj.loadDataset(obj.SMF,DatasetIndex);
            
            % Generate localizations from the current Dataset.
            LD = smi_core.LocalizeData(Dataset, obj.SMF);
            fprintf('Generating localizations ...\n');
            [SMD] = LD.genLocalizations();

            % Define NDatasets, and DatasetNum from the dataset count.
            SMD.NDatasets  = 1;
            SMD.DatasetNum = DatasetCount * ones(size(SMD.FrameNum));
            
            % Perform frame-connection on localizations in SMD.
            if obj.SMF.FrameConnection.On 
                FC = smi_core.FrameConnection(SMD, obj.SMF);
                fprintf('Frame connecting ...\n');
                [SMD] = FC.performFrameConnection();
            end
            
            % Intra-dataset drift correction.
            if obj.SMF.DriftCorrection.On 
                fprintf('Drift correcting (intra-datastet) ...\n');
                SMD = obj.DC.driftCorrectKNNIntra(SMD, DatasetIndex);
            end
        end

        % ---------------------------------------------------------------------

        function [Dataset, SMF]=loadDataset(obj,SMF,DatasetIndex)
        % loadDataset loads a dataset and converts to photons
        % set obj.Data   
        switch SMF.Data.FileType;
            case 'mat'
                [~, Dataset, SMF] = ...
                    smi_core.LoadData(SMF,SMF.Data.DataVariable,DatasetIndex);
            case 'h5'
                [~, Dataset, SMF] = smi_core.LoadData(SMF,DatasetIndex);
        end % switch
        end

        % ---------------------------------------------------------------------

        function saveResults(obj)
        % saveResults saves all results and plots in subfolder.
        % gaussblobs, drift image, fits/frame, NumConnected hist,
        % Driftcorrection plots, precision hist, intensity hist,
        % mat file with SMD and SMF files. 
        if isempty(obj.SMD)
            error('No SMD results structure found to save!');
        end

        [~, f, ~] = fileparts(obj.SMF.Data.DatasetList{1});
        if isempty(obj.SMF.Data.AnalysisID)
            fn = [f, '_Results.mat'];
        else
            fnextend = strcat('_Results_', obj.SMF.Data.AnalysisID, '.mat');
            fn = [f, fnextend];
        end

        SMD = obj.SMD;
        SMF = obj.SMF;
        save(fullfile(ResultsDir, fn), 'SMD', 'SMF', '-v7.3');
        end

    end % methods
    % =========================================================================

    % =========================================================================
    methods(Static)
        Success = unitTest();
    end % methods(Static)
    % =========================================================================
    
end
