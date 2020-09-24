classdef SMLM < handle
    % Single Molecule Localization Microscopy Analysis
    %
    % This is a high-level class that provides complete analysis of SMLM data.
    %
    %
    
    properties
        SMF
        Preset  %   {'TIRF','Sequential'} good idea?
        Data    %   Current dataset or used for manual setting of data
        DataType% {'File','UserDefined'} ?
        FileName    %String or Cell array of strings
        DataDir
        ResultsDir  % (Default = 'DataDir/../Results/FileName/) same as Seq
        
        
    end
    
    properties (Access=protected)
        SMD % SMD structure with final analysis results
        NDataSets %
        
    end
    
    
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
            
            obj.SMD=[];
            for nn=1:obj.NDatasets
                SMDnn=obj.analyzeDataset(nn);
                obj.SMD=smi_core.SingleMoleculeData.catSMD(obj.SMD,SMDnn);
            end
            
            % Inter-dataset drift correction.
            DC = smi_core.DriftCorrection(obj.SMF);
            obj.SMD = DC.driftCorrectKNNInter(obj.SMD);
        end
        
        
        function SMD=analyzeDataset(obj,DataSetIndex)
        % analyzeDataset Load and analyze one dataset    
            
            
            Dataset=obj.loadDataset(DataSetIndex);
            
            % Generate localizations from the current Dataset.
            LD = smi_core.LocalizeData(Dataset, obj.SMF);
            [SMD] = LD.genLocalizations();
            
            % Perform frame-connection on localizations in SMD.
            FC = smi_core.FrameConnection(SMD, obj.SMF);
            [SMD] = FC.performFrameConnection();
            
            % Intra-dataset drift correction.
            DC = smi_core.DriftCorrection(obj.SMF);
            SMD = DC.driftCorrectKNNIntra(SMD);
        end
        
        
        function loadDataset(obj,DataSetIndex)
        % loadDataset load a dataset and convert to photons
        % set obj.Data    
        end
        
        function saveResults(obj)
        % saveResults Save all results and plots in subfolder
        % gaussblobs, drift image, fits/frame, NumConnected hist,
        % Driftcorrection plots, precision hist, intensity hist,
        % mat file with SMD and SMF files. 
        
        
        end
        
    end
    
end
