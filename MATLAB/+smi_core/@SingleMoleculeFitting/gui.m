function gui(obj, GUIParent)
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

% Created by:
%   David J. Schodt (Lidke lab, 2020)


% Define a list of 'special' fields which are treated differently below
% (i.e., we want to manually create special GUI elements for these fields
% below).
SpecialFields = {'Data.FileName', 'Fitting.ZFitStruct'};

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
TabGroup = uitabgroup(GUIParent, 'Position', [0, 0.1, 1, 0.9], ...
    'Units', 'pixels');
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
    % initial values if appropriate.
    TopPosition = PropertyTabs{ff}.InnerPosition(2) ...
        + PropertyTabs{ff}.InnerPosition(4);
    for ss = 1:NSubfields
        % Create the 'text' and 'edit' uicontrols, skipping the 'edit'
        % control if the field is a 'struct' (an edit box wouldn't make
        % sense in that case).  If the field is a 'struct', we'll make a
        % 'popupmenu' style uicontrol for the fields.  If the field is in
        % the list of special fields defined at the top of this code, we'll
        % define it's GUI elements on a case by case basis.
        % NOTE: Some numbers were added to positions arbitrarily to 
        %       improve appearance.  
        CurrentYPosition = TopPosition - TextInitPos(4) ...
            - (UIControlInitPos(4)-0.01)*ss;
        uicontrol(PropertyTabs{ff}, 'Style', 'text', ...
            'String', [SubfieldNames{ff}{ss}, ': '], ...
            'HorizontalAlignment', 'right', 'Units', 'normalized', ...
            'Position', TextInitPos + [0, CurrentYPosition, 0, 0])
        CurrentSubfield = CurrentProperty.(SubfieldNames{ff}{ss});
        CurrentSubfieldName = [obj.SMFPropertyNames{ff}, ...
            '.', SubfieldNames{ff}{ss}];
        if ismember(CurrentSubfieldName, SpecialFields)
            % Create a listbox uicontrol to show the file names.
            if strcmp(CurrentSubfieldName, 'Data.FileName')
                % Create the listbox.
                UIControls{ff}{ss}{1} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'listbox', 'String', CurrentSubfield, ...
                    'Units', 'normalized', ...
                    'Position', UIControlInitPos ...
                    + [0, CurrentYPosition, 0, 0]);
                
                % Create a button to allow for selection of other files.
                UIControls{ff}{ss}{2} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'pushbutton', 'String', 'Select File(s)', ...
                    'Units', 'normalized', ...
                    'Position', UIControlInitPos ...
                    + [UIControlInitPos(3), CurrentYPosition, 0, 0], ...
                    'Callback', @selectFiles);
                
                % Define the location of an additional note (the note
                % contained in obj.SMFPropertyNames, to be used below).
                NotePosition = TextInitPos ...
                    + [UIControls{ff}{ss}{2}.Position(1)...
                    + UIControls{ff}{ss}{2}.Position(3), ...
                    CurrentYPosition, TextInitPos(3), 0];
            elseif strcmp(CurrentSubfieldName, 'Fitting.ZFitStruct')
                % Add an edit box associated with this pop-up menu.
                SubSubfields = fieldnames(CurrentSubfield);
                UIControls{ff}{ss}{2} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'edit', 'Units', 'normalized',...
                    'String', ...
                    num2str(CurrentSubfield.(SubSubfields{1})), ...
                    'Position', UIControlInitPos ...
                    + [UIControlInitPos(3), CurrentYPosition, 0, 0], ...
                    'Callback', @guiToProperties);
                
                % Make the pop-up menu.
                UIControls{ff}{ss}{1} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'popupmenu', ...
                    'String', fieldnames(CurrentSubfield), ...
                    'Units', 'normalized', ...
                    'Position', UIControlInitPos ...
                    + [0, CurrentYPosition, 0, 0], ...
                    'Callback', {@structCallback, ...
                    UIControls{ff}{ss}{2}, ...
                    obj.SMFPropertyNames{ff}, ...
                    SubfieldNames{ff}{ss}});
                
                % Define the location of an additional note (the note
                % contained in obj.SMFPropertyNames, to be used below).
                NotePosition = TextInitPos ...
                    + [UIControls{ff}{ss}{2}.Position(1)...
                    + UIControls{ff}{ss}{2}.Position(3), ...
                    CurrentYPosition, TextInitPos(3), 0];
            end
        else
            % Make the edit box (the default behavior for property fields).
            UIControls{ff}{ss} = uicontrol(PropertyTabs{ff}, ...
                'Style', 'edit', 'String', num2str(CurrentSubfield), ...
                'Units', 'normalized', 'Position', UIControlInitPos ...
                + [0, CurrentYPosition, 0, 0], ...
                'Callback', @guiToProperties);
            
            % Define the location of an additional note (the note
            % contained in obj.SMFPropertyNames, to be used below).
            NotePosition = TextInitPos ...
                + [UIControls{ff}{ss}.Position(1) ...
                + UIControls{ff}{ss}.Position(3), ...
                CurrentYPosition, TextInitPos(3), 0];
        end
        
        % Add the additional message for this field present in
        % obj.SMFFieldNotes.
        uicontrol(PropertyTabs{ff}, 'Style', 'text', ...
            'String', obj.SMFFieldNotes.(...
            obj.SMFPropertyNames{ff}).(SubfieldNames{ff}{ss}), ...
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
                % do something different on a case by case basis.
                CurrentFieldName = [PropertyNames{nn}, ...
                    '.', PropertyFields{mm}];
                if ismember(CurrentFieldName, SpecialFields)
                    if strcmp(CurrentFieldName, 'Data.FileName')
                        UIControls{nn}{mm}{1}.String = ...
                            ClassProperty.(PropertyFields{mm});
                    elseif strcmp(CurrentFieldName, 'Fitting.ZFitStruct')
                        UIControls{nn}{mm}{2}.String = ...
                            num2str(ClassProperty.(PropertyFields{mm}) ...
                            .(UIControls{nn}{mm}{1}.String{...
                            UIControls{nn}{mm}{1}.Value}));
                    end
                else
                    UIControls{nn}{mm}.String = ...
                        num2str(ClassProperty.(PropertyFields{mm}));
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
                    elseif strcmp(CurrentPropertyName, ...
                            'Fitting.ZFitStruct')
                        % Data.ZFitStruct has several sub-fields which we
                        % will set according to the current state of 
                        % UIControls{nn}{mm}{1}.
                        CurrentPropertyValue = obj.(PropertyNames{nn}) ...
                            .(PropertyFields{mm}) ...
                            .(UIControls{nn}{mm}{1}.String{...
                            UIControls{nn}{mm}{1}.Value});
                        obj.(PropertyNames{nn}).(PropertyFields{mm}) ...
                            .(UIControls{nn}{mm}{1}.String{...
                            UIControls{nn}{mm}{1}.Value}) = ...
                            processUserInput(...
                            UIControls{nn}{mm}{2}.String, ...
                            CurrentPropertyValue);
                    end
                else
                    % This is the default behavior for non-special fields:
                    % property is updated based on a single input in an 
                    % edit box.
                    CurrentPropertyValue = ...
                        obj.(PropertyNames{nn}).(PropertyFields{mm});
                    obj.(PropertyNames{nn}).(PropertyFields{mm}) = ...
                        processUserInput(UIControls{nn}{mm}.String, ...
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
        obj.Data.FileName = FileName.';
        
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

    function [ProcessedInput] = processUserInput(UserInput, ...
            CurrentPropertyValue)
        % This function will process a user input to an edit box uicontrol
        % so that it can be set as a property value (e.g., converting a
        % string like '12' to the number 12).
        
        % Process the input based on the existing property value.
        if isnumeric(CurrentPropertyValue)
            % NOTE: If properties should be a specific type, that should be
            %       taken care of in a set method inside 
            %       SingleMoleculeFitting.m
            % NOTE: I'm using str2num() here instead of, e.g., str2double()
            %       because str2num('') = [], but str2double('') = NaN (in
            %       MATLAB 2020b at least).
            ProcessedInput = str2num(UserInput);
        elseif islogical(CurrentPropertyValue)
            % Logical fields aren't numeric, even if they were set to
            % something like logical(1). I'm keeping this separate from
            % above condition because I anticipate the treatment of
            % numerics above to change.
            % NOTE: str2num() works nicely on logicals, e.g.,
            %       str2num('true') will be a logical (not numeric).
            ProcessedInput = str2num(UserInput);
        else
            % Note that some of the fields might be cell arrays, e.g., file
            % names.  If the user has input a string to the corresponding
            % edit box then we'll just assume they want to overwrite the
            % existing cell array.
            ProcessedInput = UserInput;
        end
    end


end