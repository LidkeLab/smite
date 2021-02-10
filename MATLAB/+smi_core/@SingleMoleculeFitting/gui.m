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
%               modified in propertiesToGUI(), guiToProperties(), and
%               perhaps in other places (look for usage of the variable
%               'SpecialFields' and the method ismember()). Note that you
%               should update and use the nested function
%               'processUserInput()' inside of propertiesToGUI() when
%               appropriate.
%           Fields that don't fit any of the categories above may or may
%               not work without change to the code below.
%       If you would like special behavior for a field, you must add that
%       field manually to the list given in a cell array in 'SpecialFields'
%       below. You must then manually add the special behavior to the
%       various appearances throughout this code (look for usage of the
%       variable 'SpecialFields' and the method ismember()).
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
    'Data.DatasetMods', 'Data.CalibrationFilePath', ...
    'Fitting.FitType', 'Fitting.ZFitStruct', 'Tracking.Method'};

% Create a figure handle for the GUI if needed.
if ~(exist('GUIParent', 'var') && ~isempty(GUIParent) ...
        && isgraphics(GUIParent))
    DefaultFigurePosition = get(0, 'defaultFigurePosition');
    GUIParent = figure('MenuBar', 'none', ...
        'Name', 'SMF Editor', 'NumberTitle', 'off', ...
        'Units', 'pixels', ...
        'Position', [DefaultFigurePosition(1:2), 485, 350]);
end

% Generate some tabs in the GUI, one per class property. While doing so,
% determine how many fields will be in each tab (we'll use this info. to
% scale certain GUI objects).
NSMFFields = numel(obj.SMFPropertyNames);
TabGroup = uitabgroup(GUIParent, 'Units', 'normalized', ...
    'Position', [0, 0.06, 1, 0.94]);
PropertyTabs = cell(NSMFFields, 1);
NFieldsMax = 0;
for ff = 1:NSMFFields
    % Create the tab.
    PropertyTabs{ff} = uitab(TabGroup, ...
        'Title', obj.SMFPropertyNames{ff}, 'units', 'normalized');
    
    % Determine how many fields will be represented in the current tab.
    NFieldsMax = max(NFieldsMax, ...
        numel(fieldnames(obj.(obj.SMFPropertyNames{ff}))));
end

% Populate each tab with useful ui features, e.g., edit boxes.
% NOTE: I've moved this into a separate loop over fields just to clean up
%       the code a bit.
TabHeight = PropertyTabs{1}.Position(4) - PropertyTabs{1}.Position(2);
TextInitPos = [0, 0, 0.25, (TabHeight/(NFieldsMax+1)) * 0.8];
UIControlInitPos = [TextInitPos(1)...
    + TextInitPos(3), 0, 0.25, TabHeight / (NFieldsMax+1)];
UIControls = cell(NSMFFields, 1);
SubfieldNames = cell(NSMFFields, 1);
for ff = 1:NSMFFields
    % Generate a list of the sub-fields to be displayed in the current tab.
    CurrentProperty = obj.(obj.SMFPropertyNames{ff});
    SubfieldNames{ff} = fieldnames(CurrentProperty);
    NSubFields = numel(SubfieldNames{ff});
    UIControls{ff} = cell(NSubFields, 1);
    
    % Create the edit boxes for simple properties (e.g., scalar
    % properties) and their associated labels, filling the edit boxes with
    % initial values if appropriate. Logical scalars will be given a
    % checkbox style uicontrol instead of an edit box. Other items
    % contained in the SpecialFields list hard-coded above will be treated
    % differently on a case-by-case basis.
    TopPosition = PropertyTabs{ff}.InnerPosition(2) ...
        + PropertyTabs{ff}.InnerPosition(4);
    for ss = 1:NSubFields
        % Create uicontrols for each sub-field. We will create either an
        % edit box or a check box (for logicals) for sub-fields unless
        % special behavior is defined (for the list of 'SpecialFields').
        % NOTE: Some numbers were added to positions arbitrarily to
        %       improve appearance.
        CurrentSubfield = CurrentProperty.(SubfieldNames{ff}{ss});
        CurrentSubfieldName = [obj.SMFPropertyNames{ff}, ...
            '.', SubfieldNames{ff}{ss}];
        CurrentFieldNote = obj.SMFFieldNotes.(...
            obj.SMFPropertyNames{ff}).(SubfieldNames{ff}{ss});
        CurrentYPosition = TopPosition - UIControlInitPos(4)*(ss+0.5);
        uicontrol(PropertyTabs{ff}, 'Style', 'text', ...
            'String', [SubfieldNames{ff}{ss}, ': '], ...
            'HorizontalAlignment', 'right', ...
            'Tooltip', CurrentFieldNote.Tip, ...
            'Units', 'normalized', ...
            'Position', TextInitPos + [0, CurrentYPosition, 0, 0])
        if ismember(CurrentSubfieldName, SpecialFields)
            % Create a listbox uicontrol to show the file names.
            switch CurrentSubfieldName
                case 'Data.FileName'
                    % Create the listbox.
                    UIControls{ff}{ss}{1} = uicontrol(PropertyTabs{ff}, ...
                        'Style', 'popupmenu', 'String', CurrentSubfield, ...
                        'Tooltip', CurrentFieldNote.Tip, ...
                        'Units', 'normalized', ...
                        'Position', UIControlInitPos ...
                        + [0, CurrentYPosition, 0, 0]);
                    
                    % Create a button to allow for selection of other
                    % files.
                    UIControls{ff}{ss}{2} = uicontrol(PropertyTabs{ff}, ...
                        'Style', 'pushbutton', ...
                        'String', 'Select File(s)', ...
                        'Tooltip', CurrentFieldNote.Tip, ...
                        'Units', 'normalized', ...
                        'Position', UIControlInitPos ...
                        + [UIControlInitPos(3), CurrentYPosition, 0, 0],...
                        'Callback', @selectRawDataFiles);
                case 'Data.CameraType'
                    % Add a pop-up menu for the CameraType options.
                    UIControls{ff}{ss} = uicontrol(PropertyTabs{ff}, ...
                        'Style', 'popupmenu', ...
                        'String', {'EMCCD', 'SCMOS'}, ...
                        'Tooltip', CurrentFieldNote.Tip, ...
                        'Units', 'normalized', ...
                        'Position', UIControlInitPos ...
                        + [0, CurrentYPosition, 0, 0], ...
                        'Callback', @guiToProperties);
                case 'Data.DatasetMods'
                    % Add a pop-up menu for the dataset modifier options.
                    UIControls{ff}{ss}{1} = uicontrol(PropertyTabs{ff}, ...
                        'Style', 'popupmenu', ...
                        'String', ...
                        {'Include Datasets', 'Exclude Datasets'}, ...
                        'Tooltip', CurrentFieldNote.Tip, ...
                        'Units', 'normalized', ...
                        'Position', UIControlInitPos ...
                        + [0, CurrentYPosition, 0, 0], ...
                        'Callback', @guiToProperties);
                    
                    % Add an edit box associated with the pop-up menu.
                    UIControls{ff}{ss}{2} = uicontrol(PropertyTabs{ff}, ...
                        'Style', 'edit', 'Units', 'normalized',...
                        'String', num2str(CurrentSubfield{1}), ...
                        'Tooltip', CurrentFieldNote.Tip, ...
                        'Position', UIControls{ff}{ss}{1}.Position ...
                        + [UIControls{ff}{ss}{1}.Position(3), 0, 0, 0], ...
                        'Callback', @guiToProperties);
                case 'Data.CalibrationFilePath'
                    % Make an edit box for the calibration file path.
                    UIControls{ff}{ss}{1} = uicontrol(PropertyTabs{ff}, ...
                        'Style', 'edit', ...
                        'String', num2str(CurrentSubfield),...
                        'Tooltip', CurrentFieldNote.Tip, ...
                        'Units', 'normalized', ...
                        'Position', UIControlInitPos ...
                        + [0, CurrentYPosition, 0, 0], ...
                        'Callback', @guiToProperties);
                    
                    % Create a button to allow for selection of the file.
                    UIControls{ff}{ss}{2} = uicontrol(PropertyTabs{ff}, ...
                        'Style', 'pushbutton', 'String', 'Select File', ...
                        'Tooltip', CurrentFieldNote.Tip, ...
                        'Units', 'normalized', ...
                        'Position',  UIControls{ff}{ss}{1}.Position ...
                        + [UIControls{ff}{ss}{1}.Position(3), 0, 0, 0], ...
                        'Callback', @selectCalDataFiles);
                case 'Fitting.FitType'
                    % Add a pop-up menu for the FitType options.
                    UIControls{ff}{ss} = uicontrol(PropertyTabs{ff}, ...
                        'Style', 'popupmenu', ...
                        'String', ...
                        {'XYNB', 'XYNBS', 'XYNBSXSY', 'XYZNB'}, ...
                        'Tooltip', CurrentFieldNote.Tip, ...
                        'Units', 'normalized', ...
                        'Position', UIControlInitPos ...
                        + [0, CurrentYPosition, 0, 0], ...
                        'Callback', @guiToProperties);
                case 'Fitting.ZFitStruct'
                    % Add an edit box associated with the pop-up menu
                    % defined below.
                    SubSubfields = fieldnames(CurrentSubfield);
                    UIControls{ff}{ss}{2} = uicontrol(PropertyTabs{ff}, ...
                        'Style', 'edit', 'Units', 'normalized',...
                        'String', ...
                        num2str(CurrentSubfield.(SubSubfields{1})), ...
                        'Tooltip', CurrentFieldNote.Tip, ...
                        'Position', UIControlInitPos ...
                        + [UIControlInitPos(3), CurrentYPosition, 0, 0],...
                        'Callback', @guiToProperties);
                    
                    % Make the pop-up menu.
                    UIControls{ff}{ss}{1} = uicontrol(PropertyTabs{ff}, ...
                        'Style', 'popupmenu', ...
                        'String', fieldnames(CurrentSubfield), ...
                        'Tooltip', CurrentFieldNote.Tip, ...
                        'Units', 'normalized', ...
                        'Position', UIControlInitPos ...
                        + [0, CurrentYPosition, 0, 0], ...
                        'Callback', @propertiesToGUI);
                case 'Tracking.Method'
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
        end
        if isempty(UIControls{ff}{ss})
            % Instead of just adding a default else condition to the
            % previous if/elseif block, we should check if
            % UIControls{ff}{ss} was even defined above. This might be
            % nice because somebody may want to add something related to
            % one of the 'SpecialFields' while still keeping the default
            % check/edit uicontrols defined here.
            if islogical(CurrentSubfield)
                % If this field is a logical type, we'll make a checkbox
                % instead of an edit box.
                UIControls{ff}{ss} = uicontrol(PropertyTabs{ff}, ...
                    'Style', 'checkbox', 'Value', CurrentSubfield, ...
                    'Tooltip', CurrentFieldNote.Tip, ...
                    'Units', 'normalized', 'Position', UIControlInitPos ...
                    + [UIControlInitPos(3)/2, CurrentYPosition, ...
                    -UIControlInitPos(3)/2, 0], ...
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
    'Tooltip', ...
    sprintf(['This button updates the GUI to reflect changes\n', ...
    'made to SMF properties outside of the GUI.']), ...
    'Units', 'normalized', 'Position', ExtraButtonsInitPos, ...
    'Callback', @propertiesToGUI);

% Add a button which can import a previously saved SMF.
uicontrol(GUIParent, 'Style', 'pushbutton', 'String', 'Import SMF', ...
    'Tooltip', ...
    sprintf(['This button allows you to import an SMF structure\n', ...
    'saved in a .mat file']), ...
    'Units', 'normalized', ...
    'Position', ExtraButtonsInitPos ...
    + [ExtraButtonsInitPos(1)+ExtraButtonsInitPos(3), 0, 0, 0],...
    'Callback', @importSMF);

% Add a button which can export the SMF (save as a struct in a .mat file).
uicontrol(GUIParent, 'Style', 'pushbutton', 'String', 'Export SMF', ...
    'Tooltip', ...
    sprintf(['This button allows you to export the current settings\n', ...
    'displayed in this GUI as a SMF structure in a .mat file']), ...
    'Units', 'normalized', ...
    'Position', ExtraButtonsInitPos ...
    + [ExtraButtonsInitPos(1) + 2*ExtraButtonsInitPos(3), 0, 0, 0],...
    'Callback', @exportSMF);

% Add a button which can reset the SMF to the default values.
uicontrol(GUIParent, 'Style', 'pushbutton', 'String', 'Reset SMF', ...
    'Tooltip', ...
    sprintf(['This button allows you to reset the current settings\n', ...
    'to their default values defined in\n', ...
    'smi_core.SingleMoleculeFitting.m']), ...
    'Units', 'normalized', ...
    'Position', ExtraButtonsInitPos ...
    + [ExtraButtonsInitPos(1) + 3*ExtraButtonsInitPos(3), 0, 0, 0],...
    'Callback', @resetSMF);

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
                CurrentField = ClassProperty.(PropertyFields{mm});
                CurrentFieldName = [PropertyNames{nn}, ...
                    '.', PropertyFields{mm}];
                if ismember(CurrentFieldName, SpecialFields)
                    switch CurrentFieldName
                        case 'Data.FileName'
                            % Data.FileName is just a cell array of
                            % strings, which is exactly what we set the
                            % uicontrol string to.
                            UIControls{nn}{mm}{1}.String = CurrentField;
                        case 'Data.DatasetMods'
                            % Data.DatasetMods has a pop-up menu defining
                            % which index of itself should be displayed
                            % (e.g., display Data.DatasetMods{2}).
                            NewString = sprintf('%i,', ...
                                CurrentField{UIControls{nn}{mm}{1}.Value});
                            UIControls{nn}{mm}{2}.String = ...
                                NewString(1:end-1);
                        case 'Data.CalibrationFilePath'
                            % Data.CalibrationFilePath should be displayed
                            % in the edit box (i.e., UIControls{nn}{mm}{1})
                            UIControls{nn}{mm}{1}.String = CurrentField;
                        case {'Data.CameraType', ...
                                'Fitting.FitType', ...
                                'Tracking.Method'}
                            % All of these fields are defined by the item
                            % selected in a pop-up menu.
                            UIControls{nn}{mm}.Value = ...
                                find(strcmp(UIControls{nn}{mm}.String, ...
                                CurrentField));
                        case 'Fitting.ZFitStruct'
                            % Fitting.ZFitStruct is defined by an edit box,
                            % but the property the edit box changes is
                            % specified by a pop-up menu.
                            UIControls{nn}{mm}{2}.String = ...
                                num2str(CurrentField ...
                                .(UIControls{nn}{mm}{1}.String{...
                                UIControls{nn}{mm}{1}.Value}));
                    end
                else
                    % Logical fields are defined by a checkbox, and we
                    % don't really care to update the 'string' of those
                    % uicontrols.
                    PropertyValue = CurrentField;
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
                    switch CurrentPropertyName
                        case 'Data.FileName'
                            % Data.FileName is just a cell array of char
                            % arrays, with the char arrays being exactly
                            % what we've set in
                            % UIControls{nn}{mm}{1}.String.
                            obj.(PropertyNames{nn}) ...
                                .(PropertyFields{mm}) = ...
                                UIControls{nn}{mm}{1}.String;
                        case 'Data.DatasetMods'
                            % Data.DatasetMods is a cell array, and we want
                            % to update the element specified by a pop-up
                            % menu.
                            CellIndex = UIControls{nn}{mm}{1}.Value;
                            CurrentPropertyValue = ...
                                obj.(PropertyNames{nn}) ...
                                .(PropertyFields{mm}){CellIndex};
                            obj.(PropertyNames{nn}) ...
                                .(PropertyFields{mm}){CellIndex} = ...
                                processUserInput(...
                                UIControls{nn}{mm}{2}, ...
                                CurrentPropertyValue);
                        case 'Data.CalibrationFilePath'
                            % Data.CalibrationFilePath can be updated
                            % directly with the string stored in
                            % UIControls{nn}{mm}{1}.
                            obj.(PropertyNames{nn}) ...
                                .(PropertyFields{mm}) = ...
                                UIControls{nn}{mm}{1}.String;
                        case {'Data.CameraType', ...
                                'Fitting.FitType', ...
                                'Tracking.Method'}
                            % All of these fields are defined by the item
                            % selected in a pop-up menu, with the string
                            % for those items being exactly the value we
                            % want to set for the class property.
                            obj.(PropertyNames{nn}) ...
                                .(PropertyFields{mm}) = ...
                                UIControls{nn}{mm}.String{ ...
                                UIControls{nn}{mm}.Value};
                        case 'Fitting.ZFitStruct'
                            % Fitting.ZFitStruct has several sub-fields
                            % which we will set according to the current
                            % state of UIControls{nn}{mm}{1}.
                            CurrentPropertyValue = ...
                                obj.(PropertyNames{nn}) ...
                                .(PropertyFields{mm}) ...
                                .(UIControls{nn}{mm}{1}.String{ ...
                                UIControls{nn}{mm}{1}.Value});
                            obj.(PropertyNames{nn})...
                                .(PropertyFields{mm}) ...
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

    function selectRawDataFiles(~, ~)
        % This is a callback function for the SMF.Data.FileName file
        % selection button of the GUI. This callback will call uiputfile()
        % to allow the user to select a set of files.
        
        % Create the file selection dialog.
        [FileName, FileDir] = uigetfile({'*.mat', 'MAT-files (*.mat)'; ...
            '*.h5', 'H5-files (*.h5)'; '*', 'Other file types'}, ...
            'Multiselect', 'on');
        if (isequal(FileName, 0) || isequal(FileDir, 0))
            return
        end
        
        % Update obj with these newly selected files.
        obj.Data.FileDir = FileDir;
        obj.Data.FileName = FileName;
        
        % Update the GUI.
        propertiesToGUI();
    end

    function selectCalDataFiles(~, ~)
        % This is a callback function for the SMF.Data.CalibrationFilePath
        % file selection button of the GUI. This callback will call
        % uiputfile() to allow the user to select the calibration file.
        
        % Create the file selection dialog.
        [FileName, FileDir] = uigetfile({'*.mat', 'MAT-files (*.mat)'; ...
            '.h5', 'H5-files (*.h5)'}, ...
            'Multiselect', 'off');
        if (isequal(FileName, 0) || isequal(FileDir, 0))
            return
        end
        
        % Update obj with these newly selected files.
        obj.Data.CalibrationFilePath = fullfile(FileDir, FileName);
        
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

    function resetSMF(~, ~)
        % This function will reset the current class instance obj to all
        % default values. This is done after first asking the user if they
        % are sure they want to proceed (accidentally clicking this would
        % be very annoying!).
        
        % Ask the user if they are sure they want to reset the SMF.
        Response = questdlg(['Are you sure you want to reset the SMF ', ...
            'properties to their defaults?'], 'Warning', ...
            'yes', 'no', 'no');
        if strcmp(Response, 'no')
            return
        end
        
        % Proceed to reset the SMF properties.
        obj.resetSMF();
        
        % Update the GUI to reflect the changes.
        propertiesToGUI();
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