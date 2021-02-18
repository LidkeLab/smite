function exportTransform(obj, FileName, FileDir)
%exportTransform exports transform information into a .mat file.
% This method will save a bunch of relevant class fields into a .mat file
% in the specified location.
%
% INPUTS:
%   FileName: Name of the .mat file we'll save. (char array/string)
%             (Default = 'RegistrationTransform' plus a string containing 
%             time)
%   FileDir: Directory in which we'll save the transform information.
%            (char array/string)(Default = obj.SMF.Data.FileDir)

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Set defaults if needed.
if (~exist('FileDir', 'var') || isempty(FileDir))
    FileDir = obj.SMF.Data.FileDir;
end
if (~exist('FileName', 'var') || isempty(FileName))
    DateString = smi_helpers.genTimeString();
    FileName = ['RegistrationTransform', '_', DateString, '.mat'];
end

% Save the registration information.
RegistrationTransform = obj.RegistrationTransform;
Coordinates = obj.Coordinates;
FiducialROI = obj.FiducialROI;
SplitFormat = obj.SplitFormat;
TransformationBasis = obj.TransformationBasis;
TransformationType = obj.TransformationType;
SeparationThreshold = obj.SeparationThreshold;
NNeighborPoints = obj.NNeighborPoints;
PolynomialDegree = obj.PolynomialDegree;
SMF = obj.SMF;
AutoscaleFiducials = obj.AutoscaleFiducials;
save(fullfile(FileDir, FileName), 'RegistrationTransform', ...
    'Coordinates',  'FiducialROI', 'SplitFormat', ...
    'TransformationBasis', 'TransformationType', ...
    'SeparationThreshold', 'NNeighborPoints', 'PolynomialDegree', ...
    'SMF', 'AutoscaleFiducials')


end