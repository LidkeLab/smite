classdef LocalizeData < handle
    %LocalizeData contains methods used to find localizations in raw data.
    % This class contains method(s) to generate localizations from numeric
    % arrays of raw data in the form of images/stacks of images.
    %
    % NOTE: All default class properties are set in the constructor.  These
    %       properties are extracted/constructed from either in constructor
    %       input SMF or from a default SMF created as
    %       SMF = smi_core.SingleMoleculeFitting.createSMF();
    %
    % EXAMPLE USAGE:
    %   LD = smi_core.LocalizeData(RawData, SMF);
    %   [SMD, SMDPreThresh] = LD.genLocalizations();
    %   OR
    %   [~, SMD, SMDPreThresh] = smi_core.LocalizeData(RawData, SMF, 1);
    %
    % REQUIRES:
    %   Image Processing Toolbox
    %   Statistics and Machine Learning Toolbox
    
    
    properties
        SMF % see smi_core.SingleMoleculeFitting
        ScaledData  % (float array)(Photons) Gain/offset corrected data
        Verbose = 1 % verbosity level
    end
    
    properties (SetAccess = 'protected')
        SMD % fully thresholded SMD, see SingleMoleculeData class
        SMDPreThresh % "full" SMD including bad localizations w/ ThreshFlag
    end
    
    methods
        function [obj, SMD, SMDPreThresh] = LocalizeData(...
                ScaledData, SMF, Verbose, AutoRun)
            %LocalizeData creates an instance of the LocalizeData class.
            % This method will prepare the LocalizeData class, setting
            % inputs to class properties if provided. This constructor can
            % also be used to run the main run method of the class, as long
            % as ScaledData and SMF are provided. To do so, set AutoRun = 1
            
            % Set a default for the AutoRun flag, which specifies whether
            % or not we should run genLocalizations() in this constructor.
            if (~exist('AutoRun', 'var') || isempty(AutoRun))
                AutoRun = 0;
            end
            
            % Set a default for the SMF structure, whose fields will be
            % extracted and set to class properties/will be used to
            % construct class properties.
            if (exist('SMF', 'var') && ~isempty(SMF))
                % Pad the input SMF structure to ensure it contains all
                % fields defined in smi_core.SingleMoleculeFitting.
                obj.SMF = smi_core.SingleMoleculeFitting.padSMF(SMF);
            else
                obj.SMF = smi_core.SingleMoleculeFitting;
            end
            
            % Set the verbosity level if provided.
            if (exist('Verbose', 'var') && ~isempty(Verbose))
                obj.Verbose = Verbose;
            end
            if ((obj.Verbose>2) && (nargin()>0))
                fprintf(['\tLocalizeData constructor: ', ...
                    'Setting class properties based on constructor ', ...
                    'inputs...\n'])
            end
                        
            % Set the ScaledData as a class property if provided.
            if (exist('ScaledData', 'var') && ~isempty(ScaledData))
                obj.ScaledData = ScaledData;
                AllFieldsSet = 1;
            else
                AllFieldsSet = 0;
            end
            
            % Run genLocalizations() if desired.
            if (AutoRun && AllFieldsSet)
                if (obj.Verbose > 2)
                    fprintf(['\tLocalizeData constructor: ', ...
                        'Auto-running obj.genLocalizations()...\n'])
                end
                [SMD, SMDPreThresh] = obj.genLocalizations();
            else
                if ((obj.Verbose>0) && (nargout()>1))
                    warning(['Constructor outputs SMD and ', ...
                        'SMDPreThresh requested but ', ...
                        'localizations were not generated.'])
                end
                SMD = [];
            end
        end
        
        [SMD, SMDPreThresh] = genLocalizations(obj);
        
    end
    
    methods (Static)
        [Success] = unitTest();
    end
    
    
end