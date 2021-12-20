function [Image] = removeBorder(Image, Border, Direction)
%removeBorder removes border pixels from the input image.
% This method isolates a central section of 'Image' which corresponds to
% deleting a border of width 'Border' from the image.
%
% INPUTS:
%   Image: Image(s) containing the border to be removed. (1-4D array)
%   Border: Border to be removed. If given as a scalar, the same border is
%           from all dimensions.
%           (scalar of NDimsx1 vector)(Default = zeros(ndims(Image), 1))
%   Direction: Border mode describing which edges are removed.
%              ('both', 'pre', 'post')(Default = 'both')
%
% OUTPUTS:
%   Image: Input 'Image' with borders removed.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults.
NDims = ndims(Image);
if (~exist('Border', 'var') || isempty(Border))
    Border = zeros(NDims, 1);
end
if (~exist('Direction', 'var') || isempty(Direction))
    Direction = 'both';
end

% Ensure that 'Border' is valid based on the size of 'Image'.
if isrow(Border)
    Border = Border.';
end
if isscalar(Border)
    % For a scalar border, the same border is removed along each dimension.
    Border = Border * ones(NDims, 1);
else
    % For a vector border, we'll assume the missing dimensions should have
    % a border of 0 (i.e., no border is removed).
    Border = padarray(Border, max(0, [NDims-numel(Border), 0]), ...
        0, 'post');
end
ImSize = size(Image);
MaxBorder = floor(ImSize / 2).';
Border = max(0, min(MaxBorder, Border(1:NDims)));

% Remove the border.
ColonIndices = repmat({':'}, NDims, 1);
for nn = 1:NDims
    NDimIndices = ColonIndices;
    switch lower(Direction)
        case 'pre'
            NDimIndices{nn} = (1+Border(nn)):ImSize(nn);
        case 'post'
            NDimIndices{nn} = 1:(ImSize(nn)-Border(nn));
        otherwise
            NDimIndices{nn} = (1+Border(nn)):(ImSize(nn)-Border(nn));
    end
    Image = Image(NDimIndices{:});
end


end