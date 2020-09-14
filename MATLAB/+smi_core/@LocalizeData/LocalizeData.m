classdef LocalizeData < handle
    %LocalizeData contains methods used to find localizations in raw data.
    
    
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
        EmitterModel % (Photons) Model of the emitter fit blobs
        ModelDataOverlay % (Photons) Overlay of ScaledData and EmitterModel
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
            
            % "Unwrap" the SMF structure to see if it contains any of the
            % class properties that we'll need.
            if (exist('SMF', 'var') && ~isempty(SMF))
                BoxFindingFields = {'BoxSize', 'MinPhotons'};
                DataFields = {'CameraType', ...
                    'CameraGain', ...
                    'CameraOffset', ...
                    'CameraReadNoise'};
                FittingFields = {'FitType', 'PSFSigma'};
                for ff = 1:numel(BoxFindingFields)
                    % Check if this field exists in SMF.BoxFinding, and if
                    % it does, copy it into the appropriate class property.
                    if isfield(SMF.BoxFinding, BoxFindingFields{ff})
                        obj.(BoxFindingFields{ff}) = ...
                            SMF.BoxFinding.(BoxFindingFields{ff});
                    end
                end
                for ff = 1:numel(DataFields)
                    % Check if this field exists in SMF.Data, and if it 
                    % does, copy it into the appropriate class property.
                    if isfield(SMF.Data, DataFields{ff})
                        obj.(DataFields{ff}) = SMF.Data.(DataFields{ff});
                    end
                end
                for ff = 1:numel(FittingFields)
                    % Check if this field exists in SMF.Fitting, and if it
                    % does, copy it into the appropriate class property.
                    if isfield(SMF.Fitting, FittingFields{ff})
                        obj.(FittingFields{ff}) = ...
                            SMF.Fitting.(FittingFields{ff});
                    end
                end
            end
        end
    end
    
    methods (Static)
        
    end
    
    
end

