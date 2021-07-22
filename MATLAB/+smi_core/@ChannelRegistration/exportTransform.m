function exportTransform(obj, TransformIndex, FileDir)
%exportTransform exports transform information into a .mat file.
% This method will save a bunch of relevant class fields into a .mat file
% in the specified location.
%
% INPUTS:
%   TransformIndex: Index of the transform that will be saved, indexed
%                   as obj.RegistrationTransform{TransformNumber}.
%   FileDir: Directory in which transforms will be saved.
%            (Default = obj.SMF.Data.FileDir)

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
    LoopIndices = TransformIndex;
else
    LoopIndices = 2:numel(obj.RegistrationTransform);
end
for nn = LoopIndices
    % Define a unique file path for this transform.
    FileName = sprintf('RegistrationTransform%ito1_%s.mat', ...
        nn, smi_helpers.genTimeString());
    
    % Save the requested transform.
    RegistrationTransform = obj.RegistrationTransform{nn};
    Coordinates = obj.Coordinates{nn};
    TransformIndex = nn;
    save(fullfile(FileDir, FileName), 'RegistrationTransform', ...
        'Coordinates',  'FiducialROI', 'SplitFormat', ...
        'TransformationBasis', 'TransformationType', ...
        'SeparationThreshold', 'NNeighborPoints', 'PolynomialDegree', ...
        'SMF', 'AutoscaleFiducials', 'TransformIndex')
end


end