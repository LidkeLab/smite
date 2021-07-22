function [FilePath] = exportTransform(obj, TransformIndex, FileDir)
%exportTransform exports transform information into a .mat file.
% This method will save a bunch of relevant class fields into a .mat file
% in the specified location.
%
% INPUTS:
%   TransformIndex: Index of the transform that will be saved, indexed
%                   as obj.RegistrationTransform{TransformNumber}.
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

% Proceed based on the provided inputs. If a 'TransformIndex' was provided, 
% we're only saving that transform. Otherwise, we'll want to save all of 
% the transforms from 2:end (transform 1 isn't useful but is kept for
% clarity in indexing).
SplitFormat = obj.SplitFormat;
TransformationBasis = obj.TransformationBasis;
TransformationType = obj.TransformationType;
SeparationThreshold = obj.SeparationThreshold;
NNeighborPoints = obj.NNeighborPoints;
PolynomialDegree = obj.PolynomialDegree;
SMF = obj.SMF;
AutoscaleFiducials = obj.AutoscaleFiducials;
FiducialROI = obj.FiducialROI;
if (obj.Verbose > 1)
    fprintf(['\tChannelRegistration.exportTransform(): ', ...
        'Exporting transform(s)...\n'])
end
if (exist('TransformIndex', 'var') && ~isempty(TransformIndex))
    TransformIndices = TransformIndex;
else
    TransformIndices = 2:numel(obj.RegistrationTransform);
end
NExports = numel(TransformIndices);
FilePath = cell(NExports, 1);
for nn = 1:NExports
    % Define a unique file path for this transform.
    FileName = sprintf('RegistrationTransform%ito1_%s.mat', ...
        TransformIndices(nn), smi_helpers.genTimeString());
    FilePath{nn} = fullfile(FileDir, FileName);
    
    % Save the requested transform.
    RegistrationTransform = obj.RegistrationTransform{TransformIndices(nn)};
    Coordinates = obj.Coordinates{TransformIndices(nn)};
    TransformIndex = TransformIndices(nn);
    save(FilePath{nn}, 'RegistrationTransform', ...
        'Coordinates',  'FiducialROI', 'SplitFormat', ...
        'TransformationBasis', 'TransformationType', ...
        'SeparationThreshold', 'NNeighborPoints', 'PolynomialDegree', ...
        'SMF', 'AutoscaleFiducials', 'TransformIndex')
end


end