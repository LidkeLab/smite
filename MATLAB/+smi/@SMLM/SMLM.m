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
            
            %interdataset drift correction
        end
        
        
        function SMD=analyzeDataset(obj,DataSetIndex)
        % analyzeDataset Load and analyze one dataset    
            
            
            Dataset=obj.loadDataset(DataSetIndex);
            
            %localizeData
            SMD=smi_core.genLocalizations(Dataset,obj.SMF);
            
            %frame connection
            
            %drift correction 
            
            
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
