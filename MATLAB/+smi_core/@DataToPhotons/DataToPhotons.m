classdef DataToPhotons < handle
    % DataToPhotons contains static methods to convert raw data to photons.
    
    
    properties
    end
    
    methods
    end
    
    methods (Static)
        [Data, ReadNoise] = convertToPhotons(RawData, SMF, ...
            RawDataROI, CalibrationROI);
    end
    
    
end

