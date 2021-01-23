classdef DataToPhotons < handle
    %DataToPhotons contains static methods to convert raw data to photons.
    %
    % This class contains static methods associated with the gain and
    % offset correction needed to convert raw data from the camera (arrays
    % given in Analog to Digital Units (ADU)) to units of photons.  The
    % main usage of this class is shown in the EXAMPLE USAGE section below.
    %
    % EXAMPLE USAGE:
    %   Given RawData in units of ADU, and an SMF structure with the fields
    %   SMF.Data.CameraGain, SMF.Data.CameraOffset, and 
    %   SMF.Data.CameraReadNoise (see SingleMoleculeFitting class for
    %   details), you can convert RawData and CameraReadNoise to units of
    %   photons and photons^2, respectively, as follows:
    %       [~, RawDataConverted, CameraReadNoiseConverted] = ...
    %           smi_core.DataToPhotons(SMF, ...
    %           RawData, RawDataROI, CalibrationROI, true);
    %   Alternatively, you can prepare the class for usage (and set class
    %   parameters immediately) as follows:
    %       DTP = smi_core.DataToPhotons(SMF, ...
    %           RawData, RawDataROI, CalibrationROI);
    
    
    properties
        % obj.RawData converted to units of photons (numeric array)
        CorrectedData
        
        % obj.CameraReadNoise converted to units of photons (numeric array)
        CorrectedReadNoise
        
        % Data that is to be gain/offset corrected (numeric array)
        RawData 
        
        % Region of interest of the raw data (numeric array)
        % (see obj.convertToPhotons() for details/usage)
        RawDataROI
        
        % Gain of the camera used to collect RawData (numeric array)
        CameraGain
        
        % Offset of the camera used to collect RawData (numeric array)
        CameraOffset
        
        % Read noise of the camera used to collect Raw Data (numeric array)
        CameraReadNoise
        
        % Region of interest of the gain/offset arrays (numeric array)
        % (see obj.convertToPhotons() for details/usage)
        CalibrationROI
    end
    
    methods
        function [obj, Data, ReadNoise] = DataToPhotons(SMF, ...
                RawData, RawDataROI, CalibrationROI, AutoRun)
            %DataToPhotons is the class constructor.
            
            % Set defaults if needed.
            if (~exist('AutoRun', 'var') || isempty(AutoRun))
                AutoRun = 0;
            end
                        
            % Set class properties based on the inputs.
            AllFieldsSet = true;
            if (exist('SMF', 'var') && ~isempty(SMF))
                % Reload the SMF to ensure it has all required properties.
                SMF = smi_core.SingleMoleculeFitting.reloadSMF(SMF);
                
                % Set the desired SMF fields to class properties.
                obj.CameraGain = SMF.Data.CameraGain;
                obj.CameraOffset = SMF.Data.CameraOffset;
                obj.CameraReadNoise = SMF.Data.CameraReadNoise;
            else
                AllFieldsSet = false;
            end
            if (exist('RawData', 'var') && ~isempty(RawData))
                obj.RawData = RawData;
            else
                AllFieldsSet = false;
            end
            if (exist('RawDataROI', 'var') && ~isempty(RawDataROI))
                obj.RawDataROI = RawDataROI;
            else
                AllFieldsSet = false;
            end
            if (exist('CalibrationROI', 'var') && ~isempty(CalibrationROI))
                obj.CalibrationROI = CalibrationROI;
            else
                AllFieldsSet = false;
            end
            
            % Run obj.convertToPhotons() if requested.
            if (AutoRun && AllFieldsSet)
                [Data, ReadNoise] = obj.convertToPhotons(obj.RawData, ...
                    obj.CameraGain, obj.CameraOffset, ...
                    obj.CameraReadNoise, ...
                    obj.RawDataROI, obj.CalibrationROI);
                obj.CorrectedData = Data;
                obj.CorrectedReadNoise = ReadNoise;
            end
        end
    end
    
    methods (Static)
        [Data, ReadNoise] = convertToPhotons(RawData, ...
            CameraGain, CameraOffset, CameraReadNoise, ...
            RawDataROI, CalibrationROI);
        [Success] = unitTest();
    end
    
    
end