function gui(obj)
%gui generates a GUI to facilitate use of the SPT class.
% This method generates a GUI which allows the user to load single-particle
% tracking data, set parameters through the SingleMoleculeFitting class,
% and generate results.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Create a figure for the GUI.
DefaultFigurePosition = get(0, 'defaultFigurePosition');
FigureXYSize = [DefaultFigurePosition(3), 600];
GUIFigure = figure('MenuBar', 'none', ...
    'Name', 'SPT Interface', 'NumberTitle', 'off', ...
    'Units', 'pixels', ...
    'Position', [DefaultFigurePosition(1), 0, FigureXYSize]);

% Add some panels to help organize the GUI.
SMFPanel = uipanel(GUIFigure, 'Title', 'Fitting parameters', ...
    'Units', 'normalized', 'Position', [0, 0.3, 1, 0.7]);
ControlPanel = uipanel(GUIFigure, 'Title', 'Controls', ...
    'Units', 'normalized', 'Position', [0, 0, 1, 0.3]);

% Stick the SingleMoleculeFitting GUI inside of the SMFPanel.
obj.SMF.gui(SMFPanel);

% Add some controls to the ControlPanel.
TextSize = [0, 0, 0.2, 0.2];
EditSize = [0, 0, 0.1, 0.2];
ButtonSize = [0, 0, 0.2, 0.2];


end