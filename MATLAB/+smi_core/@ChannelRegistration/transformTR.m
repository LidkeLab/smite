function [TRMoving] = transformTR(RegistrationTransform, TRMoving)
%transformTR transforms TR a structure using the specified transform.
% This method is just a wrapper around
% smi_core.ChannelRegistration.transformSMD().
%
% INPUTS:
%   RegistrationTransform: RegistrationTransform is a MATLAB tform object
%                          describing the transform to be visualized.
%                          (tform object)
%   TRMoving: TR structure (see smi_core.TrackingResults) that will be
%             transformed by 'RegistrationTransform'.
%
% OUTPUTS:
%   TRMoving: TR structure transformed by 'RegistrationTransform'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Make sure none of the trajectories have been transformed.  If any have,
% throw an error.
assert(~any([TRMoving.IsTransformed]), ...
    ['Some trajectories in ''TR'' has already been transformed! ', ...
    'If you intend to apply another transform to this TR, set ', ...
    'TR.IsTransformed to false for each trajectory before proceeding.']);

% Loop through trajectories in 'TRMoving' and apply the transform.
for nn = 1:numel(TRMoving)
    TRMoving(nn) = smi_core.ChannelRegistration.transformSMD(...
        RegistrationTransform, TRMoving(nn));
end


end