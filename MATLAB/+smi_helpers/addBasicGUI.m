function [ControlHandles] = addBasicGUI(GUIParent, ParamStruct, Callback)
%addBasicGUI adds simple GUI controls based on the fields in 'ParamStruct'
% This method will add basic GUI control elements to 'GUIParent', with each
% control representing a setting for a field in 'ParamStruct'.  The
% intention is that for a set of boring parameters in 'ParamStruct' (e.g.,
% scalars, single string entries, logical values, ...), you can produce a
% quick GUI to allow for setting these parameters.
%
% NOTE: This GUI does not update anything in 'ParamStruct', so all updates
%       must be done through the associated 'ControlHandles'.
%
% INPUTS:
%   GUIParent: The 'Parent' of this GUI, e.g., a figure handle.
%   ParamStruct: Structure array containing "scalar" fields (i.e., fields
%                can be set by a single edit box, checkbox, button, etc.)
%   Callback: Callback function to be set for all of the uicontrols (see
%             uicontrol property 'Callback'). (Default = '')
%             NOTE: Each of the uicontrols is given a 'Tag' property with
%                   the name of the corresponding field in 'ParamStruct',
%                   which can be useful for defining a callback function.
%
% OUTPUTS:
%   ControlHandles: Array of uicontrols corresponding to the fields in
%                   'ParamStruct', provided in the same order as the fields
%                   appear in the output of fieldnames(ParamStruct).

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Set defaults.
if ~exist('Callback', 'var')
    Callback = '';
end

% Define some parameters for the GUI controls.
FieldNames = fieldnames(ParamStruct);
NFields = numel(FieldNames);
ControlWidth = 0.5;
ControlHeight = 1 / NFields;

% Loop through the structure fields and add uicontrol elements.
ControlHandles = cell(NFields, 1);
for nn = 1:NFields
    % Isolate some info. about the nn-th field.
    CurrentName = FieldNames{nn};
    CurrentField = ParamStruct.(CurrentName);
    
    % Add a text description for the current field.
    ControlPosition = [0, 1 - nn/NFields, ControlWidth, ControlHeight];
    uicontrol('Parent', GUIParent, ...
        'Style', 'text', 'String', CurrentName, ...
        'Units', 'normalized', ... 
        'Position', ControlPosition, ...
        'HorizontalAlignment', 'left');
    
    % Create a uicontrol based on the type of the current field.
    if islogical(CurrentField)
        % Logical (boolean) fields are given a checkbox.
        ControlHandles{nn} = uicontrol(GUIParent, ...
            'Style', 'checkbox', 'Value', CurrentField, ...
            'Units', 'normalized', ...
            'Position', ControlPosition + [0.5, 0, 0, 0], ...
            'Callback', Callback, 'Tag', CurrentName);
    elseif (ischar(CurrentField) || isstring(CurrentField))
        % Chars and strings will be given an edit box.
        ControlHandles{nn} = uicontrol(GUIParent, ...
            'Style', 'edit', 'String', CurrentField, ...
            'Units', 'normalized', ...
            'Position', ControlPosition + [0.5, 0, 0, 0], ...
            'Callback', Callback, 'Tag', CurrentName);
    elseif isnumeric(CurrentField)
        % These fields are assumed to be numeric vectors, so we'll attempt 
        % to prepare edit boxes for them.
        EditString = sprintf('%g, ', CurrentField);
        ControlHandles{nn} = uicontrol(GUIParent, ...
            'Style', 'edit', 'String', EditString(1:(end-2)), ...
            'Units', 'normalized', ...
            'Position', ControlPosition + [0.5, 0, 0, 0], ...
            'Callback', Callback, 'Tag', CurrentName);
    else
        % All other fields aren't given a control since they likely can't
        % be set nicely with a simple control.
    end
end


end