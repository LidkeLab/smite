classdef LocalizeData < handle
    %LocalizeData contains methods used to find localizations in raw data.
    % This class contains method(s) to generate localizations from numeric
    % arrays of raw data in the form of images/stacks of images.
    % 
    %
    %
    % REQUIRES:
    %   DipImage, to use joinchannels() in LocalizeData.genLocalizations().
    
    properties
        CameraType = 'EMCCD'; % (Default = 'EMCCD') see GaussMLE class
        CameraGain = 1; % (Default = 1)(ADU/e-)
        CameraOffset = 0; % (Default = 0)(ADU)
        CameraReadNoise = 0; % (Default = 0)(ADU^2)
        BoxSize = 7; % (Default = 7)(Pixels) see FindROI class
        FitType = 'XYNB'; % (Default = 'XYNB') see GaussMLE class
        PSFSigma = 1.3; % (Default = 1.3)(Pixels) see GaussMLE class
        MinPhotons = 200; % (Default = 200)(Photons) see FindROI class
        RawData % (ADU) Raw data output from the camera
    end
    
    properties (SetAccess = 'protected')
        ScaledData % (Photons) Gain and offset corrected RawData
        SMD % see SingleMoleculeData class
    end
    
    methods
        function [obj] = LocalizeData(RawData, SMF)
            %LocalizeData creates an instance of the LocalizeData class.
            % This method will prepare the LocalizeData class, setting
            % inputs to class properties if provided.
            
            % Set the RawData as a class property if provided.
            if (exist('RawData', 'var') && ~isempty(RawData))
                obj.RawData = RawData;
            end
            
            % Set class properties based on the input SMF structure (if it
            % was provided).
            if (exist('SMF', 'var') && ~isempty(SMF))
                % Pad the input SMF structure to ensure it contains all
                % fields defined in SingleMoleculeFitting.createSMF().
                SMF = smi_core.SingleMoleculeFitting.padSMF(SMF);
                
                % Set the desired SMF fields to class properties.
               	obj.CameraType = SMF.Data.CameraType;
                obj.CameraGain = SMF.Data.CameraGain;
                obj.CameraOffset = SMF.Data.CameraOffset;
                obj.CameraReadNoise = SMF.Data.CameraReadNoise;
                obj.BoxSize = SMF.BoxFinding.BoxSize;
                obj.MinPhotons = SMF.BoxFinding.MinPhotons;
                obj.PSFSigma = SMF.Fitting.PSFSigma;
                obj.FitType = SMF.Fitting.FitType;
            end
        end
        
        genLocalizations(obj)
        
    end
    
    methods (Static)
        [Success] = unitTest();
    end
    
    
end

