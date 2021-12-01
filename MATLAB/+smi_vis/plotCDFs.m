function [varargout] = plotCDFs(PlotAxes, XData, CDFs, Color)
%plotCDFs plots multiple CDFs on the same plot.
% This method plots the CDFs in the cell array 'CDFs' in the same axes,
% ensuring that they both span the same range of x values (this improves
% appearance versus, e.g., having one CDF reach 1 very quickly and then not
% be present from the far right of the plot). 
%
% INPUTS:
%   PlotAxes: Axes in which CDFs will be added. (Default = gca())
%   XData: Cell array of x data corresponding to entries of CDFs.
%          (Nx1 cell array of vectors)
%   CDFs: Cell array of CDFs. (Nx1 cell array of vectors)
%   Color: Matrix defining the colors of each CDF. 
%          (Nx3 matrix of RGB triplets)(Default sampled from lines() 
%           colormap)
%
% OUTPUTS:
%   varargout{1}: Set of linehandles corresponding to the plotted CDFs.
%   varargout{2}: Axes in which the plot was made.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults and validate inputs.
if (~exist('PlotAxes', 'var') || isempty(PlotAxes) || ~isvalid(PlotAxes))
    PlotAxes = gca();
end
NCDFs = numel(CDFs);
if (~exist('Color', 'var') || isempty(Color))
    Color = lines(NCDFs);
end
assert(size(Color, 2) >= 3, 'Input ''Color'' must be a 3 column matrix!')
assert(size(Color, 1) >= NCDFs, ...
    'Input ''Color'' must have at least as many rows as the number of CDFs!')

% Pad the CDFs so that the far right tails all span the full x range
% (unless the last entry of a CDFs is < 1, in which case we'll allow it to
% stop early).  Also, make sure the data are given as column vectors.
MaxX = max(cellfun(@(X) max(X), XData));
for ii = 1:NCDFs
    % Ensure data are column vectors.
    if isrow(XData{ii})
        XData{ii} = XData{ii}.';
    end
    if isrow(CDFs{ii})
        CDFs{ii} = CDFs{ii}.';
    end
    
    % Pad the arrays.
    if ((CDFs{ii}(end)==1) && (XData{ii}(end)~=MaxX))
        XData{ii} = [XData{ii}; MaxX];
        CDFs{ii} = [CDFs{ii}; 1];
    end
end

% Plot the CDFs.
LineHandles = gobjects(NCDFs, 1);
for ii = 1:NCDFs
    LineHandles(ii) = stairs(PlotAxes, XData{ii}, CDFs{ii}, ...
        'Color', Color(ii, :));
end

% Prepare outputs.
if (nargout() > 0)
    varargout{1} = LineHandles;
end
if (nargout() > 1)
    varargout{2} = PlotAxes;
end


end