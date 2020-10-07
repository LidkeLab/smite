function gui(obj, GUIParent)
%gui is the GUI method for the SingleMoleculeFitting class.
% This method generates a GUI for the SingleMoleculeFitting class which
% allows the user to interactively view/set the class properties and the
% fields of the class properties (which are structs).
%
% NOTE: There are some properties in smi_core.SingleMoleculeFitting.m
%       whose values directly influence how this GUI works. These fields
%       should be kept up-to-date with any changes made to the
%       organization/contents of the SMF class/structures. These fields
%       include:
%           SMFPropertyNames: A cell array of char arrays/strings which
%                             will define the class properties for which we
%                             want to make a tab in this GUI.
%           SMFFieldNotes: A struct of similar organization to the SMF
%                          structures/class property organization, where
%                          the elements contain a field 'Units', which are
%                          char/string notes to be written in text boxes
%                          next to certain GUI elements. There is an
%                          additional field 'Tip', which are char/string
%                          notes to be displayed as a Tooltip when the
%                          cursor is hovered over certain GUI elements
%                          associated with that property. These notes are
%                          defined in the class constructor.
%
% NOTE: I've tried to write this to accomodate changes to the properties
%       and property sub-fields of the smi_core.SingleMoleculeFitting
%       class, but several pieces of this GUI probably still won't work
%       correctly if major changes are made. In writing this GUI this way,
%       there are a few things to keep in mind for future changes to SMF:
%           Scalar fields will automatically be given an edit box uicontrol
%           Logical fields will automatically be given a checkbox uicontrol
%           Fields whose behavior is unique should be added to the
%               hard-coded list 'SpecialFields'.  GUI elements for these
%               fields have to be manually defined below and their behavior
%               modified in propertiesToGUI, guiToProperties, and perhaps
%               in other nested functions.
%           Fields that don't fit any of the categories above may or may
%               not work without change to the code below.
%
% INPUTS:
%   GUIParent: The 'Parent' of this GUI, e.g., a figure handle.
%              (Default = figure(...))

% Created by:
%   David J. Schodt (Lidke lab, 2020)


% Define a list of 'special' fields which are treated differently below
% (i.e., we want to manually create special GUI elements for these fields
% below).
SpecialFields = {'Data.FileName', 'Data.CameraType', ...
    'Fitting.FitType', 'Fitting.ZFitStruct', 'Tracking.Method'};

% Create a figure handle for the GUI if needed.
if ~(exist('GUIParent', 'var') && ~isempty(GUIParent) ...
        && isgraphics(GUIParent))
    DefaultFigurePosition = get(0, 'defaultFigurePosition');
    GUIParent = figure('MenuBar', 'none', ...
        'Name', 'SMF Editor', 'NumberTitle', 'off', ...
        'Units', 'pixels', ...
        'Position', [DefaultFigurePosition(1:2), 490, 300]);
end

% Generate some tabs in the GUI, one per class property.
NSMFFields = numel(obj.SMFPropertyNames);
TabGroup = uitabgroup(GUIParent, 'Units', 'normalized', ...
    'Position', [0, 0.1, 1, 0.9]);
PropertyTabs = cell(NSMFFields, 1);
for ff = 1:NSMFFields
    PropertyTabs{ff} = uitab(TabGroup, ...
        'Title', obj.SMFPropertyNames{ff}, 'units', 'normalized');
end

% Populate each tab with useful ui features, e.g., edit boxes.
% NOTE: I've moved this into a separate loop over fields just to clean up
%       the code a bit.
TextInitPos = [0, 0, 0.25, 0.08];
UIControlInitPos = [TextInitPos(1)+TextInitPos(3), 0, 0.25, 0.1];
UIControls = cell(NSMFFields, 1);
SubfieldNames = UIControls;
for ff = 1:NSMFFields
    % Generate a list of the sub-fields to be displayed in the current tab.
    CurrentProperty = obj.(obj.SMFPropertyNames{ff});
    SubfieldNames{ff} = fieldnames(CurrentProperty);
    NSubfields = numel(SubfieldNames{ff});
    
    % Create the edit boxes for simple properties (e.g., scalar
    % properties) and their associated labels, filling the edit boxes with
    % initial values if appropriate. Logical scalars will be given a
    % checkbox style uicontrol instead of an edit box. Other items
    % contained in the SpecialFields list hard-coded above will be treated
    % differently on a case-by-case basis.
    TopPosition = PropertyTabs{ff}.InnerPosition(2) ...
        + PropertyTabs{ff}.InnerPosition(4);
    for ss = 1:NSubfields
        % Create the 'text' and 'edit' uicontrols, skipping the 'edit'
        % control if the field is a 'struct' (an edit box wouldn't make
        % sense in that case).  If the field is a 'struct', we'll make a
        % 'popupmenu' style uicontrol for the fields.  If the field is in
        % the list of special fields defined at the top of this code, we'll
        % define it's GUI elements on a case-by-case basis.
        % NOTE: Some numbers were added to positions arbitrarily to
        %       improve appearance.
        CurrentSubfield = CurrentProperty.(SubfieldNames{ff}{ss});
        CurrentSubfieldName = [obj.SMFPropertyNames{ff}, ...
            '.', SubfieldNames{ff}{ss}];
        CurrentFieldNote = obj.SMFFieldNotes.(...
            obj.SMFPropertyNames{ff}).(SubfieldNames{ff}{ss});
        CurrentYPosition = TopPosition - TextInitPos(4) ...
            - (UIControlInitPos(4)-0.01)*ss;
        uicontrol(PropertyTabs{ff}, 'Style', 'text', ...
            'String', [SubfieldNames{ff}{ss}, ': '], ...
            'HorizontalAlignment', 'right', ...
            'Tooltip', CurrentFieldNote.Tip, ...
            'Units', 'normalized', ...
            'Position', TextInitPos + [0, CurrentYPosition, 0, 0])
        if ismember(CurrentSubfieldName, SpecialFields)
            % Create a listbox uicontrol to show the file names.
            if strcmp(CurrentSubfieldName, 'Data.FileName')
                % Create the listbox.
                UIControls{ff}{ss}{1} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'listbox', 'String', CurrentSubfield, ...
                    'Tooltip', CurrentFieldNote.Tip, ...
                    'Units', 'normalized', ...
                    'Position', UIControlInitPos ...
                    + [0, CurrentYPosition, 0, 0]);
                
                % Create a button to allow for selection of other files.
                UIControls{ff}{ss}{2} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'pushbutton', 'String', 'Select File(s)', ...
                    'Tooltip', CurrentFieldNote.Tip, ...
                    'Units', 'normalized', ...
                    'Position', UIControlInitPos ...
                    + [UIControlInitPos(3), CurrentYPosition, 0, 0], ...
                    'Callback', @selectFiles);
            elseif strcmp(CurrentSubfieldName, 'Data.CameraType')
                % Add a pop-up menu for the CameraType options.
                UIControls{ff}{ss} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'popupmenu', ...
                    'String', {'EMCCD', 'SCMOS'}, ...
                    'Tooltip', CurrentFieldNote.Tip, ...
                    'Units', 'normalized', ...
                    'Position', UIControlInitPos ...
                    + [0, CurrentYPosition, 0, 0], ...
                    'Callback', @guiToProperties);
            elseif strcmp(CurrentSubfieldName, 'Fitting.FitType')
                % Add a pop-up menu for the FitType options.
                UIControls{ff}{ss} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'popupmenu', ...
                    'String', {'XYNB', 'XYNBS', 'XYNBSXSY', 'XYZNB'}, ...
                    'Tooltip', CurrentFieldNote.Tip, ...
                    'Units', 'normalized', ...
                    'Position', UIControlInitPos ...
                    + [0, CurrentYPosition, 0, 0], ...
                    'Callback', @guiToProperties);
            elseif strcmp(CurrentSubfieldName, 'Fitting.ZFitStruct')
                % Add an edit box associated with the pop-up menu defined
                % below.
                SubSubfields = fieldnames(CurrentSubfield);
                UIControls{ff}{ss}{2} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'edit', 'Units', 'normalized',...
                    'String', ...
                    num2str(CurrentSubfield.(SubSubfields{1})), ...
                    'Tooltip', CurrentFieldNote.Tip, ...
                    'Position', UIControlInitPos ...
                    + [UIControlInitPos(3), CurrentYPosition, 0, 0], ...
                    'Callback', @guiToProperties);
                
                % Make the pop-up menu.
                UIControls{ff}{ss}{1} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'popupmenu', ...
                    'String', fieldnames(CurrentSubfield), ...
                    'Tooltip', CurrentFieldNote.Tip, ...
                    'Units', 'normalized', ...
                    'Position', UIControlInitPos ...
                    + [0, CurrentYPosition, 0, 0], ...
                    'Callback', {@structCallback, ...
                    UIControls{ff}{ss}{2}, ...
                    obj.SMFPropertyNames{ff}, ...
                    SubfieldNames{ff}{ss}});
            elseif strcmp(CurrentSubfieldName, 'Tracking.Method')
                % Add a pop-up menu for the tracking method options.
                UIControls{ff}{ss} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'popupmenu', ...
                    'String', {'SMA_SPT'}, ...
                    'Tooltip', CurrentFieldNote.Tip, ...
                    'Units', 'normalized', ...
                    'Position', UIControlInitPos ...
                    + [0, CurrentYPosition, 0, 0], ...
                    'Callback', @guiToProperties);
            end
        else
            % If this field is a logical type, we'll make a checkbox
            % instead of an edit box.
            if islogical(CurrentSubfield)
                % Make the check box.
                UIControls{ff}{ss} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'checkbox', 'Value', CurrentSubfield, ...
                    'Tooltip', CurrentFieldNote.Tip, ...
                    'Units', 'normalized', 'Position', UIControlInitPos ...
                    + [0, CurrentYPosition, 0, 0], ...
                    'Callback', @guiToProperties);
            else
                % Make the edit box (the default behavior for property
                % fields).
                UIControls{ff}{ss} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'edit', 'String', num2str(CurrentSubfield),...
                    'Tooltip', CurrentFieldNote.Tip, ...
                    'Units', 'normalized', 'Position', UIControlInitPos ...
                    + [0, CurrentYPosition, 0, 0], ...
                    'Callback', @guiToProperties);
            end
        end
        
        % Define the location of an additional note (the note contained in
        % obj.SMFPropertyNames, to be used below).
        if iscell(UIControls{ff}{ss})
            NotePosition = TextInitPos ...
                + [UIControls{ff}{ss}{end}.Position(1)...
                + UIControls{ff}{ss}{end}.Position(3), ...
                CurrentYPosition, TextInitPos(3), 0];
        else
            NotePosition = TextInitPos ...
                + [UIControls{ff}{ss}.Position(1) ...
                + UIControls{ff}{ss}.Position(3), ...
                CurrentYPosition, TextInitPos(3), 0];
        end
        
        % Add the additional message for this field present in
        % obj.SMFFieldNotes.
        uicontrol(PropertyTabs{ff}, 'Style', 'text', ...
            'String', CurrentFieldNote.Units, ...
            'Tooltip', CurrentFieldNote.Tip, ...
            'HorizontalAlignment', 'left', 'Units', 'normalized', ...
            'Position', NotePosition)
    end
end

% Add a button which will update the GUI (could be useful if the SMF gets
% changed outside of the open GUI).
ExtraButtonsInitPos = UIControlInitPos .* [0, 0, 1, 1];
uicontrol(GUIParent, 'Style', 'pushbutton', 'String', 'Refresh GUI', ...
    'Units', 'normalized', 'Position', ExtraButtonsInitPos, ...
    'Callback', @propertiesToGUI);

% Add a button which can import a previously saved SMF.
uicontrol(GUIParent, 'Style', 'pushbutton', 'String', 'Import SMF', ...
    'Units', 'normalized', ...
    'Position', ExtraButtonsInitPos ...
    + [ExtraButtonsInitPos(1)+ExtraButtonsInitPos(3), 0, 0, 0],...
    'Callback', @importSMF);

% Add a button which can export the SMF (save as a struct in a .mat file).
uicontrol(GUIParent, 'Style', 'pushbutton', 'String', 'Export SMF', ...
    'Units', 'normalized', ...
    'Position', ExtraButtonsInitPos ...
    + [ExtraButtonsInitPos(1) + 2*ExtraButtonsInitPos(3), 0, 0, 0],...
    'Callback', @exportSMF);

% Call PropertiesToGUI here, just as a safeguard (this isn't going to do
% anything unless code has been revised above).
propertiesToGUI()

    function propertiesToGUI(~, ~)
        % This function will update the GUI with the current properties
        % present in the SingleMoleculeFitting object 'obj'.
        
        % Loop through the class properties and set the sub-fields values
        % whenever possible.
        PropertyNames = obj.SMFPropertyNames;
        for nn = 1:numel(PropertyNames)
            ClassProperty = obj.(PropertyNames{nn});
            PropertyFields = fieldnames(ClassProperty);
            for mm = 1:numel(PropertyFields)
                % Check if the current value of UIControls is itself a cell
                % array. If it is, it should have 2 elements (we assume
                % here that the first element is something like a pop-up
                % menu and the second element is the control whose string
                % we'll update). If the current field is in the
                % SpecialFields list defined at the top of this code, we'll
                % do something different on a case-by-case basis.
                CurrentFieldName = [PropertyNames{nn}, ...
                    '.', PropertyFields{mm}];
                if ismember(CurrentFieldName, SpecialFields)
                    if strcmp(CurrentFieldName, 'Data.FileName')
                        % Data.FileName is just a cell array of strings,
                        % which is exactly what we set the uicontrol string
                        % to.
                        UIControls{nn}{mm}{1}.String = ...
                            ClassProperty.(PropertyFields{mm});
                    elseif any(strcmp(CurrentFieldName, ...
                            {'Data.CameraType', 'Fitting.FitType', ...
                            'Tracking.Method'}))
                        % All of these fields are defined by the item
                        % selected in a pop-up menu.
                        UIControls{nn}{mm}.Value = ...
                            find(strcmp(UIControls{nn}{mm}.String, ...
                            ClassProperty.(PropertyFields{mm})));
                    elseif strcmp(CurrentFieldName, 'Fitting.ZFitStruct')
                        % Fitting.ZFitStruct is defined by an edit box, but
                        % the property the edit box changes is specified by
                        % a pop-up menu.
                        UIControls{nn}{mm}{2}.String = ...
                            num2str(ClassProperty.(PropertyFields{mm}) ...
                            .(UIControls{nn}{mm}{1}.String{...
                            UIControls{nn}{mm}{1}.Value}));
                    end
                else
                    % Logical fields are defined by a checkbox, and we
                    % don't really care to update the 'string' of those
                    % uicontrols.
                    PropertyValue = ClassProperty.(PropertyFields{mm});
                    if islogical(PropertyValue)
                        UIControls{nn}{mm}.Value = PropertyValue;
                    else
                        UIControls{nn}{mm}.String = num2str(PropertyValue);
                    end
                end
            end
        end
    end

    function guiToProperties(~, ~)
        % This function will update the class instance obj to reflect
        % property changes made in the GUI. This function will call
        % propertiesToGUI after making all property changes (this is done
        % both to ensure the user set value was stored as expected and to
        % show updates to "dependent" properties).
        
        % Loop through the appropriate GUI elements and update the
        % corresponding class properties.
        PropertyNames = obj.SMFPropertyNames;
        for nn = 1:numel(PropertyNames)
            PropertyFields = fieldnames(obj.(PropertyNames{nn}));
            for mm = 1:numel(PropertyFields)
                % Check if the current property is contained in the
                % SpecialFields list defined above. If it is, we'll need to
                % take special care in setting the corresponding class
                % property.
                CurrentPropertyName = [PropertyNames{nn}, ...
                    '.', PropertyFields{mm}];
                if ismember(CurrentPropertyName, SpecialFields)
                    if strcmp(CurrentPropertyName, 'Data.FileName')
                        % Data.FileName is just a cell array of char
                        % arrays, with the char arrays being exactly
                        % what we've set in
                        % UIControls{nn}{mm}{1}.String.
                        obj.(PropertyNames{nn}).(PropertyFields{mm}) = ...
                            UIControls{nn}{mm}{1}.String;
                    elseif any(strcmp(CurrentPropertyName, ...
                            {'Data.CameraType', 'Fitting.FitType', ...
                            'Tracking.Method'}))
                        % All of these fields are defined by the item
                        % selected in a pop-up menu, with the string for
                        % those items being exactly the value we want to
                        % set for the class property.
                        obj.(PropertyNames{nn}).(PropertyFields{mm}) = ...
                            UIControls{nn}{mm}.String{...
                            UIControls{nn}{mm}.Value};
                    elseif strcmp(CurrentPropertyName, ...
                            'Fitting.ZFitStruct')
                        % Fitting.ZFitStruct has several sub-fields which
                        % we will set according to the current state of
                        % UIControls{nn}{mm}{1}.
                        CurrentPropertyValue = obj.(PropertyNames{nn}) ...
                            .(PropertyFields{mm}) ...
                            .(UIControls{nn}{mm}{1}.String{...
                            UIControls{nn}{mm}{1}.Value});
                        obj.(PropertyNames{nn}).(PropertyFields{mm}) ...
                            .(UIControls{nn}{mm}{1}.String{...
                            UIControls{nn}{mm}{1}.Value}) = ...
                            processUserInput(UIControls{nn}{mm}{2}, ...
                            CurrentPropertyValue);
                    end
                else
                    % This is the default behavior for non-special fields:
                    % property is updated based on a single input in an
                    % edit box or based on the value of a checkbox (for
                    % logicals).
                    CurrentPropertyValue = ...
                        obj.(PropertyNames{nn}).(PropertyFields{mm});
                    obj.(PropertyNames{nn}).(PropertyFields{mm}) = ...
                        processUserInput(UIControls{nn}{mm}, ...
                        CurrentPropertyValue);
                end
            end
        end
        
        % Update the GUI based on the class properties.
        propertiesToGUI();
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

    function selectFiles(~, ~)
        % This is a callback function for the SMF.Data.FileName file
        % selection button of the GUI. This callback will call uiputfile()
        % to allow the user to select a set of files.
        
        % Create the file selection dialog.
        [FileName, FileDir] = uigetfile('Y:\', 'Multiselect', 'on');
        if (isequal(FileName, 0) || isequal(FileDir, 0))
            return
        end
        
        % Update obj with these newly selected files.
        obj.Data.FileDir = FileDir;
        obj.Data.FileName = FileName;
        
        % Update the GUI.
        propertiesToGUI();
    end

    function importSMF(~, ~)
        % This function will ask the user to select the desired SMF (saved
        % in a .mat file) and then call obj.importSMF() to load the
        % properties stored in that SMF.
        
        % Ask the user to select the saved SMF and load it.
        [FileName, FilePath] = uigetfile('Y:\*.mat');
        if (isequal(FileName, 0) || isequal(FilePath, 0))
            return
        end
        load(fullfile(FilePath, FileName), 'SMF')
        
        % Ensure all GUI inputs are saved to obj (the saved SMF might not
        % have all properties set, so it could make sense that the user
        % still wanted some of the manually entered GUI fields to be kept.
        guiToProperties();
        
        % Update obj with the fields present in the loaded SMF.
        obj.importSMF(SMF)
        
        % Update the GUI to reflect potential changes.
        propertiesToGUI();
    end

    function exportSMF(~, ~)
        % This function will ask the user to define the desired save
        % filename and path before calling obj.packageSMF() to create a
        % struct containing the class properties. The struct will be saved
        % as a .mat at the desired file location.
        
        % Ask the user to specify a save location.
        [FileName, FilePath] = uiputfile('Y:\SMF.mat');
        if (isequal(FileName, 0) || isequal(FilePath, 0))
            return
        end
        
        % Ensure all GUI inputs are saved to obj.
        guiToProperties();
        
        % Package the class into a struct.
        [SMF] = obj.packageSMF();
        
        % Save the SMF in the desired location.
        save(fullfile(FilePath, FileName), 'SMF')
    end

    function [ProcessedInput] = processUserInput(UIControl, ...
            CurrentPropertyValue)
        % This function will process a user input to various uicontrols and
        % processes that input as appropriate (this is meant to be a giant
        % block of if/else statements depending on several things,
        % including property type and uicontrol type).
        
        % Process the input based on the existing property value.
        if isnumeric(CurrentPropertyValue)
            % Numeric properties are expected to always be modified through
            % an edit box uicontrol.
            % NOTE: If properties should be a specific type, that should be
            %       taken care of in a set method inside
            %       SingleMoleculeFitting.m
            % NOTE: I'm using str2num() here instead of, e.g., str2double()
            %       because str2num('') = [], but str2double('') = NaN (in
            %       MATLAB 2020b at least).
            ProcessedInput = str2num(UIControl.String);
        elseif islogical(CurrentPropertyValue)
            % Logical fields are set by checkbox uicontrols, meaning that
            % we care about the 'value' of the uicontrol (not its
            % 'string').
            % NOTE: str2num() works nicely on logicals, e.g.,
            %       str2num('true') will be a logical (not numeric).
            ProcessedInput = logical(UIControl.Value);
        else
            % Note that some of the fields might be cell arrays, e.g., file
            % names.  If the user has input a string to the corresponding
            % edit box then we'll just assume they want to overwrite the
            % existing cell array.
            ProcessedInput = UIControl.String;
        end
    end


end