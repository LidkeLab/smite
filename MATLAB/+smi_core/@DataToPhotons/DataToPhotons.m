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
    %       [RawDataConverted, CameraReadNoiseConverted] = ...
    %           smi_core.DataToPhotons.convertToPhotons(RawData, SMF, ...
    %           RawDataROI, CalibrationROI);
    %   (see convertToPhotons() for detailed descriptions/usage of
    %   RawDataROI and CalibrationROI)
    
    
    properties
    end
    
    methods
    end
    
    methods (Static)
        [Data, ReadNoise] = convertToPhotons(RawData, ...
            CameraGain, CameraOffset, CameraReadNoise, ...
            RawDataROI, CalibrationROI);
        [Success] = unitTest();
    end
    
    
end