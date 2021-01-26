function [SplitROIs] = convertSplitFormatToROIs(FullROI, SplitFormat)
%convertSplitFormatToROIs converst a split format to an array of ROIs.
% This method will take a split format (e.g., obj.SplitFormat) and convert
% it to the sub-ROIs of FullROI (e.g., what we need in obj.FiducialROI).
%
% INPUTS:
%   FullROI: Full ROI of that is being split.  FullROI(3:4) must both be
%            even numbers UNLESS SplitFormat = [] or 1.
%            ([YStart, XStart, YEnd, XEnd])
%   SplitFormat: The format guideline defining how FullROI will be split
%                into sub-ROIs (integer array)(must consist of the set of
%                integers 1:numel(SplitFormat))(Default = [])
%                OPTIONS:
%                   []: This option simply returns SplitROIs = FullROI.
%                   [1]: Same action as [], but both are included for
%                        consistent usage elsewhere in the class.
%                   [1, 2]: For use when only a single fiducial image is
%                           provided.  The image will be split along
%                           columns into two equal ROIs.
%                   [1; 2]: Similar to above, but the image is split along
%                           rows.
%                   [1, 3; 2, 4]: Image is split into four equal ROIs.
%                NOTE: Entries correspond to the row indices of the output
%                      SplitROIs according to column-major indexing, e.g.,
%                      SplitROIs(3, :) corresponds to (SplitFormat==3).  
%                      As a more explicit example: 
%                      SplitFormat = [4, 3; 1, 2] means that 
%                      SplitROIs(1, :) defines the bottom left ROI
%                      SplitROIs(2, :) defines the bottom right ROI
%                      SplitROIs(3, :) defines the top right ROI
%                      SplitROIs(4, :) defines the top left ROI
%
% OUTPUTS:
%   SplitROIs: Array of ROIs corresponding to the input FullROI split as
%              guided by SplitFormat.  Each row will have the format
%              [YStart, XStart, YEnd, XEnd].  Row indices correspond to the
%              integers specifying ROIs in SplitFormat.

% Created by:
%   David J. Schodt (Lidke lab, 2021)

% Set defaults if needed.
if ~exist('SplitFormat', 'var')
    SplitFormat = [];
end

% Return the FullROI if appropriate.
NROIs = numel(SplitFormat);
if (isempty(SplitFormat) || (NROIs==1))
    SplitROIs = FullROI;
    return
end

% Revise the SplitFormat to make sure it contains all integers
% 1:numel(SplitFormat).
[~, SortIndices] = sort(SplitFormat(:));
SplitFormat(SortIndices) = 1:NROIs;

% Split the ROIs as defined by SplitFormat.
SplitROIs = NaN(NROIs, 4);
FormatSize = size(SplitFormat);
if ((FormatSize(1)==2) && (FormatSize(2)==1))
    % Split FullROI into two equal ROIs along the rows.
    SplitROIs(1, :) = [FullROI(1), FullROI(2), ...
        FullROI(3)/2, FullROI(4)];
    SplitROIs(2, :) = [FullROI(3)/2 + 1, FullROI(2), ...
        FullROI(3), FullROI(4)];
elseif ((FormatSize(1)==1) && (FormatSize(2)==2))
    % Split FullROI into two equal ROIs along the columns.
    SplitROIs(1, :) = [FullROI(1), FullROI(2), ...
        FullROI(3), FullROI(4)/2];
    SplitROIs(2, :) = [FullROI(1), FullROI(4)/2 + 1, ...
        FullROI(3), FullROI(4)];
elseif (FormatSize(1) == FormatSize(2))
    % Split FullROI into four equal ROIs.
    SplitROIs(1, :) = [FullROI(1), FullROI(2), ...
        FullROI(3)/2, FullROI(4)/2];
    SplitROIs(2, :) = [FullROI(3)/2 + 1, FullROI(2), ...
        FullROI(3), FullROI(4)/2];
    SplitROIs(3, :) = [FullROI(1), FullROI(4)/2 + 1, ...
        FullROI(3)/2, FullROI(4)];
    SplitROIs(4, :) = [FullROI(3)/2 + 1, FullROI(4)/2 + 1, ...
        FullROI(3), FullROI(4)];
else
    error(['convertSplitFormatToROIs(): size of input ', ...
        'SplitFormat does not match defined options'])
end
SplitROIs = SplitROIs(SortIndices, :);


end