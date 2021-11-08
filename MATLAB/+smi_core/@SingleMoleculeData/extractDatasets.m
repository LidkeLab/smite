function [SMD] = extractDatasets(SMD, Datasets, Compress)
%extractDatasets extracts the 'Datasets' from 'SMD'.
% This method will return the localizations in 'SMD' that arose from the
% datasets defined by integers 'Datasets'.
%
% WARNING: This method can fail if there are "constant" SMD fields that
%          happen to share the same array size as SMD.FrameNum.
%
% INPUTS:
%   SMD: A Single Molecule Data structure.
%   Datasets: Array of dataset numbers whose localizations will be
%             extracted from SMD.
%   Compress: Flag indicating dataset numbers should be compressed (e.g.,
%             datasets [2; 5; 8] will be renamed [1; 2; 3]). 
%             (Default = false)
%
% OUTPUTS:
%   SMD: A Single Molecule Data structure containing only the data from the
%        datasets defined by 'Datasets'.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Set defaults if needed.
if (~exist('Compress', 'var') || isempty(Compress))
    Compress = false;
end

% Extract the desired datasets.
SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, ...
    ismember(SMD.DatasetNum, Datasets));

% Rename datasets if needed.
if Compress
    SMD.DatasetNum = smi_helpers.compressToRange(SMD.DatasetNum);
end


end