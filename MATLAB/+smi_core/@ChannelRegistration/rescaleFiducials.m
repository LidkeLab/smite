function [ScaledFiducials] = ...
    rescaleFiducials(Fiducials, SMF, AutoScale)
%rescaleFiducials rescales the images in Fiducials as needed.
% This method "rescales" the fiducial images in Fiducials.  When
% AutoScale = false, this just means a gain/offset correction is performed.
% When AutoScale = true, the images are scaled to the range [0, 1].
%
% INPUTS:
%   Fiducials: Stack of images to be rescaled. (MxPxNImages numeric array)
%   SMF: Single Molecule Fitting structure (see
%        smi_core.SingleMoleculeFitting)
%   AutoScale: Flag to indicate whether to do gain/offset correction with
%              the gain and offset in the SMF (false), or to rescale each
%              of the NImages to the range [0, 1] (true).
%              (Boolean)(Default = false)
%
% OUTPUTS:
%   ScaledFiducials: Input images Fiducials rescaled as appropriate.
%                    (MxPxNImages numeric array)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('SMF', 'var') || isempty(SMF))
    SMF = smi_core.SingleMoleculeFitting;
end
if (~exist('AutoScale', 'var') || isempty(AutoScale))
    AutoScale = false;
end

% Rescale the fiducials as appropriate.
if AutoScale
    % Rescale the images to the range [0, 1].
    ScaledFiducials = Fiducials;
    for nn = 1:size(Fiducials, 3)
        CurrentImage = Fiducials(:, :, nn);
        ScaledFiducials(:, :, nn) = ...
            (CurrentImage-min(CurrentImage(:))) ...
            ./ max(max(CurrentImage-min(CurrentImage(:))));
    end
else
    % Perform gain/offset correction.
    [ScaledFiducials] = smi_core.DataToPhotons.convertToPhotons(...
        Fiducials, SMF.Data.CameraGain, SMF.Data.CameraOffset, []);
end


end