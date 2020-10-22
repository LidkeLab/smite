function [RegistrationTransform] = findTransform()
%findTransform finds a channel registration transform.
% This method will find a channel registration transform object that is
% intended to register coordinates from/features in the fiducial files
% specified by obj.FiducialFileNames.
%
% OUTPUTS:
%   RegistrationTransform: A MATLAB tform object describing the
%                          registration between the fiducial images.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Load the fiducial data.
LoadData = smi_core.LoadData;
NFiducials = numel(obj.FiducialFileNames);
FiducialImages = cell(NFiducials, 1);
for ii = 1:NFiducials
    [~, FiducialImages{ii}] = LoadData.loadData(obj.SMF);
end

% Proceed based on the setting of obj.TransformationBasis (which defines
% whether or not we need to find localizations from the fiducial data).
if strcmp(obj.TransformationBasis, 'coords')
    % Use smi_core.LocalizeData to find localizations in the fiducial
    % files.
else
end




end