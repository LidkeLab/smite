function [FigHandle] = plotCumDrift(SMD, FieldName)
%plotCumDrift creates cumulative plots for Drift in any direction
% plotCumDrift will create a cumulative drift plot for the specified
% direction given the data in the SMD (single molecule results) structure
% and a CumName specifying the drift field within SMD which will be
% plotted.
%
% INPUT:
%    SMD: Single molecule results structure containing that results from
%         a previous single molecule analysis.
%    FieldName: The field name of the drift within SMD, e.g. 'DriftX'.  
%               FieldName can be either a character array, a string, or a 
%               cell array containing multiple character arrays or strings 
%               (to plot multiple drift lines on the same plot).
%
% OUTPUT:
%    FigHandle:  Figure handle of cumulative Drift plot.

% Hanieh Mazloom-Farsibaf, April 2017 (Keith A. Lidke's lab)
% David James Schodt, December 2018 (Keith A. Lidke's lab)

% Check inputs and set defaults if needed.
if ~exist('SMD', 'var')
    error('You must enter an SMD structure to plot the cumulative drift.')
end
if ~exist('FieldName', 'var')
    error('You must enter the SMD field name for the drift of interest.')
end

% Extract the drift vector(s) from the SMD structure.
if iscell(FieldName)
    % Initialize a cell array for the drift arrays.
    Drift = cell(numel(FieldName), 1);
    
    % Populate the drift cell array, ensuring each element is a vector by
    % stacking the vectors for each dataset.
    for ii = 1:numel(FieldName)
        Drift{ii} = SMD.(FieldName{ii})(:);
    end
else
    % Only one FieldName was given, so we just need one drift array (which
    % we'll place in a cell array for consistency).
    Drift = {SMD.(FieldName)(:)};
end

% Plot the drift array(s) specified by FieldName.
FigHandle = figure;
hold('on')
for ii = 1:numel(Drift)
    % Loop through each of the drift vectors specified in FieldName.
    Frames = (1:numel(Drift{ii})).';
    plot(Frames, Drift{ii}, '.')
end

% Add labels to the figure (titles, axis labels, etc.).
% NOTE: I'm turning off legend warnings in case one of SMD.(FieldName{ii})
%       is empty.  This shouldn't cause errors elsewhere so it's easy to
%       just turn off the warning temporarily.
title('Cumulative Drift')
xlabel('Absolute Frame Number')
ylabel('Shift (Camera Pixels)')
warning('off', 'MATLAB:legend:IgnoringExtraEntries')
legend(FieldName)
warning('on', 'MATLAB:legend:IgnoringExtraEntries')

end
