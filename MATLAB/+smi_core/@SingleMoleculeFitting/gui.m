function [GUIParent] = gui(obj, GUIParent)
%gui is the GUI method for the SingleMoleculeFitting class.
% This method generates a GUI for the SingleMoleculeFitting class which
% allows the user to interactively view/set the class properties and the
% fields of the class properties (which are structs).
%
% NOTE: I've tried to write this to accomodate changes to the properties
%       and property sub-fields of the smi_core.SingleMoleculeFitting 
%       class, but several pieces of this GUI probably still won't work
%       correctly if major changes are made.  For example, the selected
%       figure size of the default GUIParent might not be big enough to
%       accomodate the addition of large numbers of sub-fields to one of
%       the properties (the result would be GUI controls outside of the
%       viewable part of the figure!).
%
% INPUTS:
%   GUIParent: The 'Parent' of this GUI, e.g., a figure handle.
%              (Default = figure(...))
%   
% OUTPUTS:
%   GUIParent: The 'Parent' of this GUI, which will either be a figure
%              handle generated in this method or the input Parent provided
%              by the user.

% Created by:
%   David J. Schodt (Lidke lab, 2020)


% Create a figure handle for the GUI if needed. 
if ~(exist('GUIParent', 'var') && ~isempty(GUIParent) ...
        && isgraphics(GUIParent))
    DefaultFigurePosition = get(0, 'defaultFigurePosition');
    GUIParent = figure('MenuBar', 'none', ...
        'Name', 'SMF Editor', 'NumberTitle', 'off', ...
        'Units', 'pixels', ...
        'Position', [DefaultFigurePosition(1:2), 500, 300]);
end

% Generate some tabs in the GUI, one per class property.
NSMFFields = numel(obj.SMFFieldNames);
TabGroup = uitabgroup(GUIParent, 'Position', [0, 0.05, 1, 0.95], ...
    'Units', 'pixels');
PropertyTabs = cell(NSMFFields, 1);
for ff = 1:NSMFFields
    PropertyTabs{ff} = uitab(TabGroup, 'Title', obj.SMFFieldNames{ff}, ...
        'units', 'pixels');
end

% Populate each tab with useful ui features, e.g., edit boxes.
% NOTE: I've moved this into a separate loop over fields just to clean up
%       the code a bit.
UIControlInitPos = [100, 0, 100, 20];
TextInitPos = [0, 0, 100, 25];
UIControls = cell(NSMFFields, 1);
SubfieldNames = UIControls;
for ff = 1:NSMFFields
    % Generate a list of the sub-fields to be displayed in the current tab.
    CurrentProperty = obj.(obj.SMFFieldNames{ff});
    SubfieldNames{ff} = fieldnames(CurrentProperty);
    NSubfields = numel(SubfieldNames{ff});
    
    % Create the edit boxes for simple properties (e.g., scalar
    % properties) and their associated labels, filling the edit boxes with
    % initial values if appropriate.
    TopPosition = PropertyTabs{ff}.InnerPosition(2) ...
        + PropertyTabs{ff}.InnerPosition(4);
    for ss = 1:NSubfields
        % Create the 'text' and 'edit' uicontrols, skipping the 'edit'
        % control if the field is a 'struct' (an edit box wouldn't make
        % sense in that case).  If the field is a 'struct', we'll make a
        % 'popupmenu' style uicontrol for the fields. Note that the
        % hard-coded numbers (e.g., multiplicative factor of 2), were
        % chosen arbitrarily to improve appearance.
        CurrentYPosition = TopPosition - UIControlInitPos(4)*ss ...
            - 2*TextInitPos(4);
        uicontrol(PropertyTabs{ff}, 'Style', 'text', ...
            'String', sprintf('%s:', SubfieldNames{ff}{ss}), ...
            'HorizontalAlignment', 'right', ...
            'Position', TextInitPos + [0, CurrentYPosition-5, 0, 0])
        CurrentSubfield = CurrentProperty.(SubfieldNames{ff}{ss});
        if isstruct(CurrentSubfield)
            % Add an edit box associated with this pop-up menu.
            SubSubfields = fieldnames(CurrentSubfield);
            UIControls{ff}{ss}{2} = uicontrol(PropertyTabs{ff}, ...
                'Style', 'edit', ...
                'String', num2str(CurrentSubfield.(SubSubfields{1})), ...
                'Position', UIControlInitPos ...
                + [UIControlInitPos(1), CurrentYPosition, 0, 0]);
            
            % Make the pop-up menu.
            UIControls{ff}{ss}{1} = uicontrol(PropertyTabs{ff}, ...
                'Style', 'popupmenu', ...
                'String', fieldnames(CurrentSubfield), ...
                'Position', UIControlInitPos ...
                + [0, CurrentYPosition, 0, 0], ...
                'Callback', {@structCallback, ...
                UIControls{ff}{ss}{2}, ...
                obj.SMFFieldNames{ff}, ...
                SubfieldNames{ff}{ss}});
        else
            % Make the edit box.
            UIControls{ff}{ss} = uicontrol(PropertyTabs{ff}, ...
                'Style', 'edit', 'String', num2str(CurrentSubfield),...
                'Position', UIControlInitPos ...
                + [0, CurrentYPosition, 0, 0]);
        end
    end
end

% Add a button which will update the GUI (could be useful if the SMF gets
% changed outside of the open GUI).
uicontrol(GUIParent, 'Style', 'pushbutton', 'String', 'Refresh GUI', ...
    'Position', UIControlInitPos + [-UIControlInitPos(1:2), 0, 0], ...
    'Callback', {@propertiesToGUI, UIControls});

% Call PropertiesToGUI here, just as a safeguard (this isn't going to do
% anything unless code has been revised above).
propertiesToGUI([], [], UIControls)

    function propertiesToGUI(~, ~, UIControls)
        % This function will update the GUI with the current properties
        % present in the SingleMoleculeFitting object 'obj'.
        
        % Loop through the class properties and set the sub-fields values
        % whenever possible.
        PropertyNames = obj.SMFFieldNames;
        for nn = 1:numel(PropertyNames)
            ClassProperty = obj.(PropertyNames{nn});
            PropertyFields = fieldnames(ClassProperty);
            for mm = 1:numel(PropertyFields)
                % Check if this field is a structure. If it is, we don't
                % want to update it here (the intention is that structures
                % should have other uicontrols which have callbacks, and
                % that the GUI is updated in the callback).
                if ~isstruct(ClassProperty.(PropertyFields{mm}))
                    UIControls{nn}{mm}.String = ...
                        num2str(ClassProperty.(PropertyFields{mm}));
                end
            end
        end
    end

    function updateClassProperty(Source, ~, PropertyName, FieldName)
        % This is a callback for several uicontrols which will pass along a
        % user input contained in Source.String.  This user input will be
        % stored in obj.(PropertyName).(FieldName), converting type if
        % needed.
        
        % Set the specified class property as appropriate (based on the
        % 'type' of that property).
        CurrentPropertyValue = obj.(PropertyName).(FieldName);
        if isnumeric(CurrentPropertyValue)
        else
            % Note that some of the fields might be cell arrays, e.g., file
            % names.  If the user has input a string to the corresponding
            % edit box 
        end
    end

    function structCallback(Source, ~, ...
            AssociatedUIControl, PropertyName, FieldName)
        % This function will update the GUI anytime special struct related
        % uicontrols are interacted with (e.g., for the pop-up menus,
        % clicking one of the pop-up items should call this method). The
        % input AssociatedUIControl will be something like an edit box
        % whose 'String' property will be updated based on the 'Value'
        % property of the calling uicontrol.
        AssociatedUIControl.String = num2str(...
            obj.(PropertyName).(FieldName).(Source.String{Source.Value}));
    end


end