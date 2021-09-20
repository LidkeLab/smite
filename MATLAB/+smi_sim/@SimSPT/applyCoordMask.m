function [Coordinates, MaskedCoordinates] = ...
    applyCoordMask(Coordinates, Mask, FrameSize)
%applyCoordMask applies a (discrete) mask to the provided coordinates.
% This method applies the binary mask in 'Mask' to the input 'Coordinates',
% meaning that all coordinates falling outside of the mask are thrown out.
%
% INPUTS:
%   Coordinates: XY coordinates to be masked. ([Y, X])
%   Mask: Logical array where coordinates within true pixels are kept, and
%         coordinates within false pixels are thrown out. 
%         (YSizexXSize logical array)
%   FrameSize: Size of the frame corresponding to the coordinates in
%              'Coordinates'. This can be used in conjunction with a 'Mask'
%              which has smaller pixels than the coordinates (e.g., 'Mask'
%              could be a 128x128 matrix but coordinates came from
%              FrameSize = [32, 32], allowing for finer scale masking
%              of the coordinates). (Default = size(Mask))
%
% OUTPUTS:
%   Coordinates: Input 'Coordinates' after the 'Mask' has been applied.
%   MaskedCoordinates: Entries of input 'Coordinates' that were masked to
%                      produce the output 'Coordinates'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults and check inputs.
if (~exist('Mask', 'var') || isempty(Mask))
    return
end
CoordSize = size(Coordinates);
if ((CoordSize(2)>CoordSize(1)) && (CoordSize(1)==2))
    % This looks like an array organized as [Y; X] instead of the desired
    % [Y, X].
    Coordinates = Coordinates.';
end
MaskSize = size(Mask);
if (~exist('FrameSize', 'var') || isempty(FrameSize))
    FrameSize = MaskSize;
end

% Convert the coordinates to a binary locations on a grid matching the
% input 'Mask'.
BinaryLocations = round(Coordinates .* (MaskSize./FrameSize));
BinaryLocations = max(1, min(BinaryLocations, MaskSize));
LocationIndices = sub2ind(MaskSize, ...
    BinaryLocations(:, 1), BinaryLocations(:, 2));
WithinMask = logical(Mask(LocationIndices));
Coordinates = Coordinates(WithinMask, :);
MaskedCoordinates = Coordinates(~WithinMask, :);


end