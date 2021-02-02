function [SMDMoving, SMDFixed] = transformSMD(...
            RegistrationTransform, SMDMoving, SMDFixed)
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
%   SMDFixed: SMD structure (see SingleMoleculeData class) for the 
%             "fixed" or "reference" fiducial (see
%             ChannelRegistration.transformCoords()).
%
% OUTPUTS:
%   SMDMoving: SMD structure (see SingleMoleculeData class) for the 
%              "moving" fiducial (see
%              ChannelRegistration.transformCoords()) which can be
%              considered transformed after passing through this method
%              (the actual transformed SMD will be either SMDMoving or 
%              SMDFixed, depending on the transform used).
%   SMDFixed: SMD structure (see SingleMoleculeData class) for the 
%             "SMDFixed" fiducial (see
%             ChannelRegistration.transformCoords()) which can be
%             considered transformed after passing through this method
%             (the actual transformed SMD will be either SMDMoving or 
%             SMDFixed, depending on the transform used).

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Call transformCoords on the input SMDs.
MovingCoordinates = [SMDMoving.X, SMDMoving.Y];
FixedCoordinates = [SMDFixed.X, SMDFixed.Y];
MovingCoordsTransformed = ...
    smi_core.ChannelRegistration.transformCoords(...
    RegistrationTransform, MovingCoordinates);
SMDMoving.X = MovingCoordsTransformed(:, 1);
SMDMoving.Y = MovingCoordsTransformed(:, 2);
SMDFixed.X = FixedCoordinates(:, 1);
SMDFixed.Y = FixedCoordinates(:, 2);


end