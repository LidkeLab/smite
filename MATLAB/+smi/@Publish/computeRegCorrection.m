function [RegCorrection] = computeRegCorrection(SMF)
%computeRegCorrection computes the registration corrections made.
% This method loads the channel registration results in the raw data file
% specified in 'SMF.Data' and computes the total registration correction
% made per dataset.
%
% INPUTS:
%   FilePath: Full filepath to the raw data .h5 file.
%
% OUTPUTS:
%   RegCorrection: Registration correction made per dataset.
%                  NOTE: This is saved in units of micrometers in the raw
%                        data file.  It is converted to pixels using the
%                        value stored in SMF.Data.PixelSize.
%                  (Pixels)(NDatasetsx1)

% Created by:
%   David J. Schodt (Lidke Lab 2021)


% Load the registration results.
AlignReg = smi_core.LoadData.readH5File(...
    fullfile(SMF.Data.FileDir, SMF.Data.FileName{1}), ...
    'AlignReg');

% Compute the total correction made before each dataset.
RegCorrection = NaN(numel(AlignReg), 3);
RegCorrection(1, :) = sum(AlignReg(1).Data.ErrorSignalHistory, 1);
for ii = 2:numel(AlignReg)
    RegCorrection(ii, :) = sum(AlignReg(ii).Data.ErrorSignalHistory, 1) ...
        - sum(RegCorrection(1:(ii-1), :), 1);
end
RegCorrection = RegCorrection / SMF.Data.PixelSize;


end