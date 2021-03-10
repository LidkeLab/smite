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
    
    
    properties
        CameraType  % (char array) see GaussMLE class for usage
        BoxSize     % (scalar)(Pixels) see FindROI class for usage
        BoxOverlap  % (scalar)(Pixels) see FindROI class for usage
        FitType     % (char array) see GaussMLE class for usage/options
        PSFSigma    % (scalar)(Pixels) see GaussMLE class
        MinPhotons  % (scalar)(Photons) see FindROI class
        Iterations  % (scalar) see GaussMLE class
        ZFitStruct  % (struct array) see GaussMLE class
        MinMax      % (struct array), see Threshold class for usage
        ScaledData  % (float array)(Photons) Gain/offset corrected data
        Verbose = 1 % verbosity level
        ThresholdingOn % (true/false) perform thresholding
    end
    
    properties (SetAccess = 'protected')
        SMD % fully thresholded SMD, see SingleMoleculeData class
        SMDPreThresh % "full" SMD including bad localizations w/ ThreshFlag
    end
    
    methods
        function [obj, SMD, SMDPreThresh] = LocalizeData(...
                ScaledData, SMF, AutoRun)
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
                SMF = smi_core.SingleMoleculeFitting.padSMF(SMF);
            else
                SMF = smi_core.SingleMoleculeFitting;
            end
            
            % Set class properties based on the SMF structure input to the
            % constructor/set to a default above.
            obj.CameraType = SMF.Data.CameraType;
            obj.BoxSize = SMF.BoxFinding.BoxSize;
            obj.BoxOverlap = SMF.BoxFinding.BoxOverlap;
            obj.MinPhotons = SMF.BoxFinding.MinPhotons;
            obj.FitType = SMF.Fitting.FitType;
            obj.PSFSigma = SMF.Fitting.PSFSigma;
            obj.Iterations = SMF.Fitting.Iterations;
            obj.ZFitStruct = SMF.Fitting.ZFitStruct;
            obj.ThresholdingOn = SMF.Thresholding.On;
            obj.MinMax.X_SE = [0, SMF.Thresholding.MaxXY_SE];
            obj.MinMax.Y_SE = [0, SMF.Thresholding.MaxXY_SE];
            obj.MinMax.Z_SE = [0, SMF.Thresholding.MaxZ_SE];
            obj.MinMax.PValue = [SMF.Thresholding.MinPValue, 1];
            obj.MinMax.PSFSigma = [SMF.Thresholding.MinPSFSigma, ...
                SMF.Thresholding.MaxPSFSigma];
            obj.MinMax.Photons = [SMF.Thresholding.MinPhotons, inf];
            obj.MinMax.Bg = [0, SMF.Thresholding.MaxBg];
            
            % Set the ScaledData as a class property if provided.
            if (exist('ScaledData', 'var') && ~isempty(ScaledData))
                obj.ScaledData = ScaledData;
                AllFieldsSet = 1;
            else
                AllFieldsSet = 0;
            end
            
            % Run genLocalizations() if desired.
            if (AutoRun && AllFieldsSet)
                [SMD, SMDPreThresh] = obj.genLocalizations();
            else
                if (nargout() > 1)
                    warning(['Constructor outputs SMD and ', ...
                        'SMDPreThresh requested but ', ...
                        'localizations were not generated.'])
                end
                SMD = [];
            end
        end
        
        [SMD, SMDPreThresh] = genLocalizations(obj)
        
    end
    
    methods (Static)
        [Success] = unitTest();
    end
    
    
end
