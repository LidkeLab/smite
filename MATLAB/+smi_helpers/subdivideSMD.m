function [SMDSub, SMDROIs] = subdivideSMD(SMD, SubROISize)
%subdivideSMD divvys up an SMD into sub-ROIs.
% This method divvys up the input localizations in 'SMD' into sub-ROIs 
% matching the size specified in SubROISize, with edge cases taken to be 
% the largest ROI possible (that is smaller than the nominal sub-ROI size).
%
% NOTE: The way I've divided the ROI allows for some localizations to be
%       present in multiple sub-ROIs.
%
% INPUT:
%   SMD: Single Molecule Data structure containing the localizations
%        that'll be subdivided.
%   SubROISize: The nominal size of the output subdivided SMDs.
%               (Pixels)(2x1 array)(Default = [SMD.YSize, SMD.XSize])
%
% OUTPUT:
%   SMDSub: Structure of the subdivided SMDs.
%   SMDROIs: ROIs of the regions corresponding to the divided SMDs.
%            (NROIsx4 array)([YStart, XStart, YEnd, XEnd])
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
DataSize = [SMD.YSize, SMD.XSize];
if (~exist('SubROISize', 'var') || isempty(SubROISize))
    SubROISize = DataSize.';
end

% Ensure arrays are properly shaped.
if isrow(SubROISize)
    SubROISize = SubROISize.';
end

% Define the ROIs of the subdivided images.
NDivisions = ceil(DataSize.' ./ SubROISize);
YStart = repmat(1 + SubROISize(1)*(0:(NDivisions(1)-1)).', ...
    [NDivisions(2), 1]);
XStart = repelem(1 + SubROISize(2)*(0:(NDivisions(2)-1)).', ...
    NDivisions(1));
SMDROIs = [YStart, XStart, ...
    min(DataSize(1), YStart+SubROISize(1)-1), ...
    min(DataSize(2), XStart+SubROISize(2)-1)];

% Loop through each ROI and isolate that section of the image.
NROIs = prod(NDivisions);
SMDSub = SMD;
for nn = NROIs:-1:1
    SMDSub(nn, 1) = smi_core.SingleMoleculeData.isolateSubSMD(SMD, ...
        (ceil(SMD.Y)>=SMDROIs(nn, 1)) & (ceil(SMD.X)>=SMDROIs(nn, 2)) ...
        & (floor(SMD.Y)<=SMDROIs(nn, 3)) & (floor(SMD.X)<=SMDROIs(nn, 4)));
end


end