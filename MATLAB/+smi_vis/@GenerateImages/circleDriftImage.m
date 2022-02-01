function [CircleDriftImage, SRImageZoom] = ...
    circleDriftImage(SMD, SRImageZoom, ...
    MinPixelsPerCircle, SEScaleFactor)
%circleImageDrift generates a drift image with circles for localizations.
% This method generates a circle drift image of the localizations in SMR.
% At each localization in SMD, a circle will be placed in the image
% centered at the localization coordinates, with the radius being the mean
% of X_SE and Y_SE for that localization.
%
% INPUTS:
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData).
%   SRImageZoom: Zoom factor of the output images w.r.t. the coordinate
%                system of the localizations.
%                (scalar, integer)(Default = 20)
%   MinPixelsPerCircle: Approximately the number of pixels used for the
%                       smallest SE circle. Note that SRImageZoom takes
%                       precedence over this parameter, so you must set
%                       SRImageZoom = [] for this to work.
%                       (Default = 16 for guidance, but isn't used!)
%   SEScaleFactor: Multiplicative scaling of the mean X/Y standard errors.
%                  This is useful if you have very low SEs but still want
%                  their circles visible, without the need to increase
%                  SRImageZoom. (Default = 1);
%
% OUTPUTS:
%   CircleDriftImage: Circle image with colors representing time, as
%                     defined by the parula() colormap.
%                     (single array)(dynamic range of [0, 1])
%   SRImageZoom: Same as input 'SRImageZoom' unless the requested zoom
%                produces an image too large to be written with imwrite().
%                In that case, this will be the value which was actually
%                used.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Define default parameters.
if (~exist('SRImageZoom', 'var') ...
        || (isempty(SRImageZoom)&&isempty(MinPixelsPerCircle)))
    % NOTE: The special condition here is to ensure SRImageZoom takes
    %       precedence over MinPixelsPerCircle when appropriate (see
    %       documentation for INPUTS: MinPixelsPerCircle above).
    SRImageZoom = 20;
end
if (~exist('MinPixelsPerCircle', 'var') || isempty(MinPixelsPerCircle))
    % As it's currently written, this default won't actually be used (see
    % relation to 'SRImageZoom' and default setting.
    MinPixelsPerCircle = 16;
end
if (~exist('SEScaleFactor', 'var') || isempty(SEScaleFactor))
    SEScaleFactor = 1;
end

% Prepare the drift image.
ColorMap = parula(SMD.NFrames * SMD.NDatasets);
[~, CircleDriftImage] = smi_vis.GenerateImages.circleImage(...
    SMD, ColorMap(SMD.NFrames*(SMD.DatasetNum-1) + SMD.FrameNum, :), ...
    SRImageZoom, MinPixelsPerCircle, SEScaleFactor);


end