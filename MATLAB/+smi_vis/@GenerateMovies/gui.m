function gui(obj)
%gui prepares a movie generation GUI for the GenerateMovies class.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Create a figure handle for the GUI.
DefaultFigurePosition = get(0, 'defaultFigurePosition');
obj.GUIFigure = figure('MenuBar', 'none', ...
    'Name', 'Movie maker', 'NumberTitle', 'off', ...
    'Units', 'pixels', ...
    'Position', [DefaultFigurePosition(1:2), 1000, 500]);

% Add a panel for display options.
ParamsPanelPos = [0, 0, 0.25, 1];
ParamsPanel = uipanel(obj.GUIFigure, ...
    'Title', 'Parameters', ...
    'FontWeight', 'bold', 'Units', 'normalized', ...
    'Position', ParamsPanelPos);
smi_helpers.addBasicGUI(ParamsPanel, obj.Params, @updateParams);

% Add some controls for the movie.
ControlPanelPos = [ParamsPanelPos(3)-ParamsPanelPos(1), 0, 0.5, 0.1];
ControlPanel = uipanel(obj.GUIFigure, 'Title', 'Movie Controls', ...
    'FontWeight', 'bold', 'Units', 'normalized', ...
    'Position', ControlPanelPos);
FrameSliderPos = [0.1, 0.5, 0.89, 0.5];
FrameSlider = uicontrol('Parent', ControlPanel, ...
    'Style', 'slider', ...
    'Units', 'normalized', 'Position', FrameSliderPos, ...
    'HorizontalAlignment', 'left', ...
    'Min', 1, 'Max', obj.TRInternal(1).NFrames, ...
    'SliderStep', [1, 1] / obj.TRInternal(1).NFrames, ...
    'Value', 1, ...
    'Callback', @frameSlider);
FrameNumDisplayPos = [FrameSliderPos(1), 0, ...
    FrameSliderPos(3), 1 - FrameSliderPos(4)];
FrameNumDisplay = uicontrol('Parent', ControlPanel, ...
    'Style', 'text', ...
    'Units', 'normalized', 'Position', FrameNumDisplayPos, ...
    'HorizontalAlignment', 'left');
PlayButtonPos = [0, FrameSliderPos(2), ...
    0.9 * (1-FrameSliderPos(3)), FrameSliderPos(4)];
uicontrol('Parent', ControlPanel, ...
    'Style', 'pushbutton', 'String', 'PLAY', ...
    'Units', 'normalized', 'Position', PlayButtonPos, ...
    'Callback', @playMovieCallback);

% Prepare axes for the movie.
MoviePanelPos = ControlPanelPos ...
    + [0, ControlPanelPos(4), 0, 1-2*ControlPanelPos(4)];
MoviePanel = uipanel(obj.GUIFigure, ...
    'Units', 'normalized', 'Position', MoviePanelPos);
obj.MovieAxes = axes(MoviePanel);
obj.MovieAxes.ActivePositionProperty = 'position';
axtoolbar(obj.MovieAxes, 'default');

% Add controls to allow for saving a movie.
SaveMoviePanelPos = [MoviePanelPos(1)+MoviePanelPos(3), 0, ...
    1-ParamsPanelPos(3)-MoviePanelPos(3), ControlPanelPos(4)];
SaveMoviePanel = uipanel(obj.GUIFigure, 'Title', 'Save Movie', ...
    'FontWeight', 'bold', 'Units', 'normalized', ...
    'Position', SaveMoviePanelPos);
uicontrol('Parent', SaveMoviePanel, 'Style', 'pushbutton', ...
    'String', 'Save Movie', ...
    'Units', 'normalized', 'Position', [0, 0, 1, 1], ...
    'Callback', @saveMovieButtonClicked);

% Add a panel to contain trajectory information about clicked trajectories.
CurrentTrajIndex = min([obj.TRInternal.ConnectID]);
TrajInfoPanelPos = ...
    [SaveMoviePanelPos(1), sum(SaveMoviePanelPos([2, 4])), ...
    SaveMoviePanelPos(3), MoviePanelPos(4)];
TrajInfoPanel = uipanel(obj.GUIFigure, ...
    'Title', 'Trajectory Information', ...
    'FontWeight', 'bold', 'Units', 'normalized', ...
    'Position', TrajInfoPanelPos);
TrajInfoTitlePos = [0, 0.9, 1, 0.1];
uicontrol('Parent', TrajInfoPanel, 'Style', 'text', ...
    'Units', 'normalized', 'Position', TrajInfoTitlePos, ...
    'String', ...
    sprintf('(click a trajectory to display trajectory information)'));
TrajInfoTextPos = [0, TrajInfoTitlePos(2)-TrajInfoTitlePos(4), 0.5, 0.05];
uicontrol('Parent', TrajInfoPanel, 'Style', 'text', ...
    'String', 'Connect ID: ', ...
    'Units', 'normalized', ...
    'Position', TrajInfoTextPos, ...
    'HorizontalAlignment', 'left');
TrajIDDisplay = uicontrol('Parent', TrajInfoPanel, 'Style', 'edit', ...
    'Units', 'normalized', ...
    'Position', TrajInfoTextPos + [TrajInfoTextPos(3), 0, 0, 0], ...
    'HorizontalAlignment', 'left', ...
    'Callback', @trajectorySelectedEdit);
StartFrameDisplay = uicontrol('Parent', TrajInfoPanel, 'Style', 'text', ...
    'String', sprintf('Start %s (%s):', ...
    obj.TimeDimensionString, obj.TimeUnitString), ...
    'Units', 'normalized', ...
    'Position', TrajInfoTextPos ...
    + [0, -TrajInfoTextPos(4), (1-TrajInfoTextPos(3)), 0], ...
    'HorizontalAlignment', 'left');
EndFrameDisplay = uicontrol('Parent', TrajInfoPanel, 'Style', 'text', ...
    'String', sprintf('End %s (%s):', ...
    obj.TimeDimensionString, obj.TimeUnitString), ...
    'Units', 'normalized', 'Position', TrajInfoTextPos ...
    + [0, -2*TrajInfoTextPos(4), (1-TrajInfoTextPos(3)), 0], ...
    'HorizontalAlignment', 'left');
uicontrol('Parent', TrajInfoPanel, 'Style', 'pushbutton', ...
    'String', 'Display Plots', ...
    'Units', 'normalized', 'Position', [0, 0, 1, 0.2], ...
    'Callback', @displayPlots);

% If raw data is available, display the first frame.
if ~isempty(obj.RawData)
    obj.prepRawData()
    obj.prepAxes()
    obj.makeFrame(obj.MovieAxes, ...
        obj.TRInternal, obj.ScaledData(:, :, :, end), ...
        obj.Params, obj.SMF, obj.SMD, obj.Params.ZFrames(1));
end


    function updateParams(Source, ~)
        % Callback for obj.Params updates.
        
        % Update the appropriate field of obj.Params based on the
        % uicontrol.
        CurrentField = obj.Params.(Source.Tag);
        if islogical(CurrentField)
            obj.Params.(Source.Tag) = Source.Value;
        elseif (ischar(CurrentField) || isstring(CurrentField))
            obj.Params.(Source.Tag) = Source.String;
        elseif isnumeric(CurrentField)
            % NOTE: We want to use str2num() instead of, e.g., str2double()
            %       because str2num() works nicely for vectors.
            obj.Params.(Source.Tag) = str2num(Source.String);
        end
        
        % Update the frame slider to reflect changes to 'ZFrames'.
        if ~isempty(obj.Params.ZFrames)
            FrameSlider.Value = obj.Params.ZFrames(1);
            FrameSlider.Min = obj.Params.ZFrames(1);
            FrameSlider.Max = obj.Params.ZFrames(2);
        end
    end

    function playMovieCallback(Source, ~)
        % This is a callback function for the PLAY button, which will call
        % the function playMovie to replay the movie without needing to use
        % the slidebar.
        
        try
            % Disable the PLAY button so it can't be clicked multiple times
            Source.Enable = 'off';
            
            % Play the movie.
            obj.generateMovie();
            
            % Remake the last frame of the movie showing the entirety of
            % the trajectories.
            TempParams = obj.Params;
            TempParams.MaxTrajLength = inf;
            obj.LineHandles = obj.makeFrame(obj.MovieAxes, ...
                obj.TRInternal, obj.ScaledData(:, :, :, end), TempParams, ...
                obj.SMF, obj.SMD, TempParams.ZFrames(2));
            
            % Set the callback function for the line handles.
            setLineCallbacks()
        catch MException
            % Re-enable the Play button before throwing the error.
            Source.Enable = 'on';
            rethrow(MException)
        end
    end

    function frameSlider(Source, ~)
        % This is a callback function for the event that the user has slid
        % the FrameSlider slidebar.  This method will determine the
        % location of the slidebar, convert this to a frame within the
        % movie, and then plot the trajectories of the movie as they exists
        % up through that frame.

        % Get the location of the slide bar, as well as the min and max
        % values possible for the slider.
        SliderValue = round(get(Source, 'Value'));
        SliderMax = get(Source, 'Max');
        
        % Update the scaled data based on the current obj.Params.
        obj.prepRawData()
        
        % Make sure the axes are prepared based on the settings in
        % obj.Params (this can be a bit slow since each movement of the
        % slider calls this, but I haven't found a better way to deal with
        % this since updates in obj.setVitalParams() can affect the changes
        % needed to the axes).
        obj.prepAxes()
        
        % Display the selected movie frame.  If the selected frame is the
        % last frame, display the entire trajectory.
        if (SliderValue == SliderMax)
            TempParams = obj.Params;
            TempParams.MaxTrajLength = inf;
            obj.LineHandles = obj.makeFrame(obj.MovieAxes, obj.TRInternal, ...
                obj.ScaledData(:, :, :, SliderValue), ...
                TempParams, obj.SMF, obj.SMD, SliderValue);
        else
            obj.LineHandles = obj.makeFrame(obj.MovieAxes, obj.TRInternal, ...
                obj.ScaledData(:, :, :, SliderValue), ...
                obj.Params, obj.SMF, obj.SMD, SliderValue);
        end
        
        % Set the callback function for the line handles.
        setLineCallbacks()

        % Update the frame number display within the GUI to show the
        % currently selected frame of the movie.
        FrameNumDisplay.String = sprintf(...
            'Frame %i of %i', SliderValue, SliderMax);
    end

    function setLineCallbacks()
        % Loop through the line handles and set their 'Callback' property
        % so that they can be clicked.
        for nn = 1:numel(obj.LineHandles)
            % Specify the ButtonDownFcn.
            if isgraphics(obj.LineHandles(nn)) ...
                    && obj.LineHandles(nn).isvalid
                % If the nn-th line handle has been deleted/is not
                % valid, we can't set the ButtonDownFcn property.
                obj.LineHandles(nn).ButtonDownFcn = ...
                    {@trajectoryClicked, nn};
            end
        end
    end

    function trajectoryClicked(~, ~, TRIndex)
        % This is a callback function for the event that the user has
        % clicked on a trajectory (with TR index 'TRIndex') within
        % the movie.
        
        % Remove previous indicator(s) that a trajectory had been selected
        % prior to this event, ensuring that the line is still valid before
        % attempting to modify it.
        [obj.LineHandles(obj.LineHandles.isvalid ...
            & isgraphics(obj.LineHandles)).LineWidth] = deal(0.5);

        % Add indicators to the selected trajectory to emphasize which one
        % was selected.
        if (CurrentTrajIndex == TRIndex)
            % If the same trajectory was clicked again, we'll assume the
            % user was trying to "unclick" it.
            CurrentTrajIndex = [];
            return
        end
        [obj.LineHandles(TRIndex).LineWidth] = 3;

        % Update 'CurrentlySelectedTrajectory' so that other callbacks can 
        % access this information.
        CurrentTrajIndex = TRIndex;

        % Grab useful information about the trajectory to be displayed,
        % ensuring the units are consistent with the movie.
        FrameNum = obj.TRInternal(TRIndex).FrameNum;
        Time = (FrameNum-1)*obj.Params.UnitFlag/obj.SMF.Data.FrameRate ...
            + FrameNum*(~obj.Params.UnitFlag);

        % Modify text panels within the display panel to provide
        % information about the selected trajectory.
        TrajIDDisplay.String = num2str(obj.TRInternal(TRIndex).ConnectID);
        StartFrameDisplay.String = sprintf('Start %s (%s): %i', ...
            obj.TimeDimensionString, ...
            obj.TimeUnitString, min(Time));
        EndFrameDisplay.String = sprintf('End %s (%s): %i', ...
            obj.TimeDimensionString, ...
            obj.TimeUnitString, max(Time));
    end

    function trajectorySelectedEdit(Source, ~)
        % This is a callback function for the event that a user has
        % manually entered a trajectory ID that they wish to highlight/get
        % information for.

        % Ensure that the entered value isn't empty.
        if isempty(Source.String)
            return
        end

        % Determine the trajectory ID and check if a line handle exists for
        % it.
        ConnectID = str2double(Source.String);
        if ~ismember(ConnectID, [obj.TRInternal.ConnectID])
            warning('Trajectory %i does not exist in the movie', ...
                ConnectID);
            return
        end

        % Call the trajectorySelected function, passing along the user
        % input trajectory ID.
        trajectoryClicked([], [], find([obj.TRInternal.ConnectID] == ConnectID))
    end

    function saveMovieButtonClicked(~, ~)
        % This is a callback function to respond to clicks of the save
        % movie button.
        obj.saveMovie()
    end


    function displayPlots(~, ~)
        % This is a callback function to respond to clicks of the Display
        % Plots button.  Clicking that button will open a small GUI with a
        % dropdown menu to allow for display of various plots related to
        % the currently clicked trajectory.
        
        % Create the figure for the small GUI.
        DisplayGUI = figure('NumberTitle', 'off', 'Resize', 'off', ...
            'Units', 'pixels', 'MenuBar', 'none', ...
            'ToolBar', 'none', 'Position', [500, 500, 200, 150]);
        
        % Add text to the GUI to guide user control.
        uicontrol('Parent', DisplayGUI, 'Style', 'text', ...
            'String', 'Select Plot:', 'Position', [70, 125, 60, 15]);
        
        % Add a drop-down menu to the GUI figure.
        DropdownMenu = uicontrol('Parent', DisplayGUI, ...
            'Style', 'popupmenu', 'String', obj.DispPlotsOptions, ...
            'Position', [50, 95, 100, 25]);
        
        % Add a display button to the GUI figure.
        uicontrol('Parent', DisplayGUI, 'Style', 'pushbutton', ...
            'String', 'Display Plot', 'Position', [65, 45, 70, 25], ...
            'Callback', {@displayPlotPushed, DropdownMenu});
    end


    function displayPlotPushed(~, ~, DropdownMenu)
        % This is a callback for the Display Plot button within the display
        % plot GUI, which will plot the information about a clicked
        % trajectory from the drop down menu.

        % Create a new figure and plot the trajectory information as
        % selected by the drop down menu in the display plot GUI.
        DisplayPlotsFigure = findobj('Tag', 'DisplayPlotsFigure');
        if isempty(DisplayPlotsFigure)
            DisplayPlotsFigure = figure('Tag', 'DisplayPlotsFigure');
        end
        figure(DisplayPlotsFigure)
        PlotAxes = axes(DisplayPlotsFigure);
        hold(PlotAxes, 'on');
        
        % Extract the desired field from the TR structure.
        DropdownString = DropdownMenu.String;
        DropdownValue = DropdownMenu.Value;
        DesiredPlotString = DropdownString{DropdownValue};
        DesiredField = obj.TRInternal(CurrentTrajIndex).(DesiredPlotString);
        FrameNum = obj.TRInternal(CurrentTrajIndex).FrameNum;
        assert(numel(FrameNum) == numel(DesiredField), ...
            ['The selected field is not the same length as ', ...
            'TR.FrameNum and cannot be plotted'])
        plot(PlotAxes, FrameNum, DesiredField, 'x')
        xlabel(PlotAxes, 'Frame number (frames)')
    end


end