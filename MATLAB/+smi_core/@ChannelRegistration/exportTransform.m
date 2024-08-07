function [FilePath] = exportTransform(obj, FileDir)
%exportTransform exports transform information into a .mat file.
% This method will save a bunch of relevant class fields into a .mat file
% in the specified location.
%
% INPUTS:
%   FileDir: Directory in which transforms will be saved.
%            (Default = obj.SMF.Data.FileDir)
%
% OUTPUTS:
%   FilePath: Cell array of path(s) to the exported transforms.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Set defaults if needed.
if (~exist('FileDir', 'var') || isempty(FileDir))
    FileDir = obj.SMF.Data.FileDir;
end
if (exist('TransformIndex', 'var') && ~isempty(TransformIndex))
    TransformIndices = TransformIndex;
else
    TransformIndices = 1:numel(obj.RegistrationTransform);
end

% Compute the registration errors (useful to export here, even if it slows
% things down).
RegistrationError = zeros(numel(TransformIndices), 1);
RegistrationErrorLOO = zeros(numel(TransformIndices), 1);
RegParams = {smi_helpers.arrayMUX({obj.PolynomialDegree, obj.NNeighborPoints}, ...
    strcmpi(obj.TransformationType, 'lwm'))};
for ii = TransformIndices
    RegistrationError(ii) = sqrt(mean(obj.estimateRegistrationError(...
        obj.RegistrationTransform{ii}, ...
        obj.Coordinates{ii}(:, :, 2), obj.Coordinates{ii}(:, :, 1))));
    RegistrationErrorLOO(ii) = sqrt(mean(obj.estimateRegErrorLOO(...
        obj.TransformationType, RegParams, ...
        obj.Coordinates{ii}(:, :, 2), obj.Coordinates{ii}(:, :, 1))));
end

% Save the requested transform(s).
if (obj.Verbose > 1)
    fprintf(['\tChannelRegistration.exportTransform(): ', ...
        'Exporting transform(s)...\n'])
end
SplitFormat = obj.SplitFormat;
TransformationBasis = obj.TransformationBasis;
TransformationType = obj.TransformationType;
SeparationThreshold = obj.SeparationThreshold;
NNeighborPoints = obj.NNeighborPoints;
PolynomialDegree = obj.PolynomialDegree;
SMF = obj.SMF;
AutoscaleFiducials = obj.AutoscaleFiducials;
FiducialROI = obj.FiducialROI;
[~, FiducialFileName] = fileparts(obj.SMF.Data.FileName{1});
ExportFileName = ...
    sprintf('RegistrationTransform_%s.mat', FiducialFileName);
FilePath = fullfile(FileDir, ExportFileName);
RegistrationTransform = obj.RegistrationTransform(TransformIndices);
Coordinates = obj.Coordinates{TransformIndices};
save(FilePath, 'RegistrationTransform', ...
    'Coordinates',  'FiducialROI', 'SplitFormat', ...
    'TransformationBasis', 'TransformationType', ...
    'SeparationThreshold', 'NNeighborPoints', 'PolynomialDegree', ...
    'SMF', 'AutoscaleFiducials', 'TransformIndices', ...
    'RegistrationError', 'RegistrationErrorLOO')


end