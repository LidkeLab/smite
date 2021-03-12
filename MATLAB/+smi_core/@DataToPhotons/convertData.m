function [CorrectedData, CorrectedReadNoise] = convertData(obj)
%convertData performs gain/offset correction on data.
% This method is meant to be a wrapper around convertToPhotons(), providing
% a means to do the gain/offset correction while storing all parameters,
% raw data, and corrected data in the class instance obj.
%
% OUTPUTS:
%   CorrectedData: gain/offset corrected obj.RawData. (float array)
%   CorrectedReadNoise: gain/offset corrected obj.ReadNoise

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Call obj.convertToPhotons() and store the outputs in obj.
if (obj.Verbose > 0)
    fprintf(['\tDataToPhotons.convertData(): ', ...
        'Performing gain and offset correction on provided data...\n'])
end
[CorrectedData, CorrectedReadNoise] = ...
    obj.convertToPhotons(obj.RawData, ...
    obj.CameraGain, obj.CameraOffset, obj.CameraReadNoise, ...
    obj.RawDataROI, obj.CalibrationROI);
obj.CorrectedData = CorrectedData;
obj.CorrectedReadNoise = CorrectedReadNoise;
if (obj.Verbose > 0)
    fprintf(['\tDataToPhotons.convertData(): ', ...
        'Gain and offset correction complete.\n'])
end


end