function [SMDMoving] = transformSMD(RegistrationTransform, SMDMoving)
%transformSMD transforms SMD structures using the specified transform.
% This method is just a wrapper around
% smi_core.ChannelRegistration.transformCoords(...) which was added for
% convenience.
%
% INPUTS:
%   RegistrationTransform: RegistrationTransform is a MATLAB tform object
%                          describing the transform to be visualized.
%                          (tform object)
%   SMDMoving: SMD structure (see SingleMoleculeData class) for the
%              "moving" fiducial (see
%              ChannelRegistration.transformCoords()).
%
% OUTPUTS:
%   SMDMoving: SMD structure (see SingleMoleculeData class) for the
%              "moving" fiducial (see
%              ChannelRegistration.transformCoords()) which can be
%              considered transformed after passing through this method
%              (the actual transformed SMD will be either SMDMoving or
%              SMDFixed, depending on the transform used).

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Check if the provided SMD has already been transformed.  If it has, don't
% proceed.
if SMDMoving.IsTransformed
    error(['The input ''SMD'' has already been transformed! If you intend ', ...
        'to apply another transform to this SMD, set SMD.IsTransformed ', ...
        'to false before proceeding.']);
end

% Call transformCoords on the input SMDs.
MovingCoordinates = [SMDMoving.X, SMDMoving.Y];
MovingCoordsTransformed = ...
    smi_core.ChannelRegistration.transformCoords(...
    RegistrationTransform, MovingCoordinates);
SMDMoving.X = MovingCoordsTransformed(:, 1);
SMDMoving.Y = MovingCoordsTransformed(:, 2);

% Update the IsTransformed field.
SMDMoving.IsTransformed = true;


end