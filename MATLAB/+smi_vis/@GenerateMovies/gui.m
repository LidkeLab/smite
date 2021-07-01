function gui(obj)
%gui prepares a movie generation GUI for the GenerateMovies class.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Create a figure handle for the GUI.
DefaultFigurePosition = get(0, 'defaultFigurePosition');
PlotFigure = figure('MenuBar', 'none', ...
    'Name', 'Movie maker', 'NumberTitle', 'off', ...
    'Units', 'pixels', ...
    'Position', [DefaultFigurePosition(1:2), 1000, 500]);

% Add a panel for display options.
ParamsPanelPos = [0, 0, 0.25, 1];
ParamsPanel = uipanel(PlotFigure, ...
    'Title', 'Parameters', ...
    'FontWeight', 'bold', 'Units', 'normalized', ...
    'Position', ParamsPanelPos);
ParamsControls = smi_helpers.addBasicGUI(ParamsPanel, obj.Params, ...
    @updateParams);

% Add some controls for the movie.
ControlPanelPos = [ParamsPanelPos(3), 0, 0.5, 0.1];
ControlPanel = uipanel(PlotFigure, 'Title', 'Movie Controls', ...
    'FontWeight', 'bold', 'Units', 'normalized', ...
    'Position', ControlPanelPos);
FrameSliderPos = [0.1, 0.5, 0.89, 0.5];
FrameSlider = uicontrol('Parent', ControlPanel, ...
    'Style', 'slider', ...
    'Units', 'normalized', 'Position', FrameSliderPos, ...
    'HorizontalAlignment', 'left', ...
    'Min', 1, 'Max', obj.TR(1).NFrames, ...
    'SliderStep', [1, 1]/obj.TR(1).NFrames, ...
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
PlayButton = uicontrol('Parent', ControlPanel, ...
    'Style', 'pushbutton', 'String', 'PLAY', ...
    'Units', 'normalized', 'Position', PlayButtonPos, ...
    'Callback', @playMovieCallback);

% Prepare axes for the movie.
MoviePanelPos = ControlPanelPos ...
    + [0, ControlPanelPos(4), 0, 1-2*ControlPanelPos(4)];
MoviePanel = uipanel(PlotFigure, ...
    'Units', 'normalized', 'Position', MoviePanelPos);
obj.MovieAxes = axes(MoviePanel);
obj.MovieAxes.ActivePositionProperty = 'position';
axtoolbar(obj.MovieAxes, 'default');

% Add controls to allow for saving a movie.
SaveMoviePanelPos = [MoviePanelPos(1)+MoviePanelPos(3), 0, ...
    1-ParamsPanelPos(3)-MoviePanelPos(3), ControlPanelPos(4)];
SaveMoviePanel = uipanel(PlotFigure, 'Title', 'Save Movie', ...
    'FontWeight', 'bold', 'Units', 'normalized', ...
    'Position', SaveMoviePanelPos);
uicontrol('Parent', SaveMoviePanel, 'Style', 'pushbutton', ...
    'String', 'Save Movie', ...
    'Units', 'normalized', 'Position', [0, 0, 1, 1], ...
    'Callback', @saveMovieButtonClicked);

% Add a panel to contain trajectory information about clicked trajectories.
CurrentlySelectedTrajectory = [];
TrajInfoPanelPos = ...
    [SaveMoviePanelPos(1), sum(SaveMoviePanelPos([2, 4])), ...
    SaveMoviePanelPos(3), MoviePanelPos(4)];
TrajInfoPanel = uipanel(PlotFigure, ...
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
    'String', 'Trajectory ID: ', ...
    'Units', 'normalized', ...
    'Position', TrajInfoTextPos, ...
    'HorizontalAlignment', 'left');
uicontrol('Parent', TrajInfoPanel, 'Style', 'edit', ...
    'Units', 'normalized', ...
    'Position', TrajInfoTextPos + [TrajInfoTextPos(3), 0, 0, 0], ...
    'HorizontalAlignment', 'left', ...
    'Callback', @trajectorySelectedEdit, 'Tag', 'TrajIDDisplay');
uicontrol('Parent', TrajInfoPanel, 'Style', 'text', ...
    'String', sprintf('Start %s (%s):', ...
    obj.TimeDimensionString, obj.TimeUnitString), ...
    'Units', 'normalized', ...
    'Position', TrajInfoTextPos ...
    + [0, -TrajInfoTextPos(4), (1-TrajInfoTextPos(3)), 0], ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'StartFrameDisplay');
uicontrol('Parent', TrajInfoPanel, 'Style', 'text', ...
    'String', sprintf('End %s (%s):', ...
    obj.TimeDimensionString, obj.TimeUnitString), ...
    'Units', 'normalized', 'Position', TrajInfoTextPos ...
    + [0, -2*TrajInfoTextPos(4), (1-TrajInfoTextPos(3)), 0], ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'EndFrameDisplay');
uicontrol('Parent', TrajInfoPanel, 'Style', 'pushbutton', ...
    'String', 'Display Plots', ...
    'Units', 'normalized', 'Position', [0, 0, 1, 0.2], ...
    'Callback', @displayPlots);

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
    end

    function playMovieCallback(Source, ~)
        % This is a callback function for the PLAY button, which will call
        % the function playMovie to replay the movie without needing to use
        % the slidebar.

%         % Disable the PLAY button so it can't be clicked multiple times.
%         Source.Enable = 'off';

        % Clear the movie axes to ensure old plots aren't still remaining.
        cla(obj.MovieAxes);

        % Play the movie.
        obj.generateMovie();
        
        % Remake the last frame of the movie showing the entirety of the
        % trajectories.
        TempParams = obj.Params;
        TempParams.MaxTrajLength = inf;
        obj.LineHandles = obj.makeFrame(obj.MovieAxes, ...
            obj.TR, obj.ScaledData(:, :, end), TempParams, obj.SMD, ...
            TempParams.ZFrames(2));
        
        % Set the callback function for the line handles.
        setLineCallbacks()
        
%         % Re-enable the PLAY button.
%         Source.Enable = 'on';
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
        
        % Ensure some needed properties are populated and updated based on
        % the current obj.Params.
        obj.rescaleData()
        
        % Make sure the axes are prepared based on the settings in
        % obj.Params (this can be a bit slow since each movement of the
        % slider calls this, but I haven't found a better way to deal with
        % this since updates in obj.setVitalParams() can affect the changes
        % needed to the axes).
        obj.prepAxes()
        
        % Display the selected movie frame.
        obj.makeFrame(obj.MovieAxes, obj.TR, ...
            obj.ScaledData(:, :, SliderValue), ...
            obj.Params, obj.SMD, SliderValue)

        % Update the frame number display within the GUI to show the
        % currently selected frame of the movie.
        FrameNumDisplay.String = sprintf(...
            'Frame %i of %i', SliderValue, SliderMax);
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
        if (CurrentlySelectedTrajectory == TRIndex)
            % If the same trajectory was clicked again, we'll assume the
            % user was trying to "unclick" it.
            CurrentlySelectedTrajectory = [];
            return
        end
        [obj.LineHandles(TRIndex).LineWidth] = 3;

        % Update 'CurrentlySelectedTrajectory'  so that other callbacks can 
        % access this information.
        CurrentlySelectedTrajectory = TRIndex;

        % Grab useful information about the trajectory to be displayed,
        % ensuring the units are consistent with the movie.
        FrameNum = obj.TR(TRIndex).FrameNum;
        Time = (FrameNum-1)*obj.Params.UnitFlag/obj.SMF.Data.FrameRate ...
            + FrameNum*(~obj.Params.UnitFlag);
% 
%         % Modify text panels within the display panel to provide
%         % information about the selected trajectory.
%         for kk = 1:numel(TrajInfoStruct.TrajInfoPanel.Children)
%             switch TrajInfoStruct.TrajInfoPanel.Children(kk).Tag
%                 case 'TrajIDDisplay'
%                     TrajInfoStruct.TrajInfoPanel.Children(kk).String = ...
%                         num2str(ii);
%                 case 'StartFrameDisplay'
%                     TrajInfoStruct.TrajInfoPanel.Children(kk).String = ...
%                         sprintf('Start %s (%s): %i', ...
%                         TrajInfoStruct.TimeDimensionString, ...
%                         TrajInfoStruct.TimeUnitString, min(Time));
%                 case 'EndFrameDisplay'
%                     TrajInfoStruct.TrajInfoPanel.Children(kk).String = ...
%                         sprintf('End %s (%s): %i', ...
%                         TrajInfoStruct.TimeDimensionString, ...
%                         TrajInfoStruct.TimeUnitString, max(Time));
%             end
%         end
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
                    {@trajectoryClicked, obj.TR(nn).ConnectID};
            end
        end
    end
%
%
%     function trajectorySelectedEdit(Source, ~)
%         % This is a callback function for the event that a user has
%         % manually entered a trajectory ID that they wish to highlight/get
%         % information for.
%
%         % Ensure that the entered value isn't empty.
%         if isempty(Source.String)
%             return;
%         end
%
%         % Determine the trajectory ID and check if a line handle exists for
%         % it.
%         TrajectoryID = str2double(Source.String);
%         if ~ismember(TrajectoryID, TrajInfoStruct.LineHandleIDMap)
%             warning('Trajectory %i does not exist in the movie', ...
%                 TrajectoryID);
%             return;
%         end
%
%         % Call the trajectorySelected function, passing along the user
%         % input trajectory ID.
%         trajectorySelected(TrajectoryID)
%     end
%

%
%     function displayPlots(~, ~)
%         % This is a callback function to respond to clicks of the Display
%         % Plots button.  Clicking that button will open a small GUI with a
%         % dropdown menu to allow for display of various plots related to
%         % the currently clicked trajectory.
%
%         % Create the figure for the small GUI.
%         DisplayGUI = figure('NumberTitle', 'off', 'Resize', 'off', ...
%             'Units', 'pixels', 'MenuBar', 'none', ...
%             'ToolBar', 'none', 'Position',[500, 500, 200, 150]);
%
%         % Add text to the GUI to guide user control.
%         uicontrol('Parent', DisplayGUI, 'Style', 'text', ...
%             'String', 'Select Plot:', 'Position', [70, 125, 60, 15]);
%
%         % Add a drop-down menu to the GUI figure.
%         PlotDropdownMenu = uicontrol('Parent', DisplayGUI, ...
%             'Style', 'popupmenu', 'String',  {'Photons'}, ...
%             'Position', [50, 95, 100, 25]);
%
%         % Add a display button to the GUI figure.
%         uicontrol('Parent', DisplayGUI, 'Style', 'pushbutton', ...
%             'String', 'Display Plot', 'Position', [65, 45, 70, 25], ...
%             'Callback', @displayPlotPushed);
%
%         % Store the handle to the dropdown menu in the TrajInfoStruct (this
%         % might be a sloppy way to pass this around but I don't see it
%         % causing any issues).
%         TrajInfoStruct.PlotDropdownMenu = PlotDropdownMenu;
%     end
%
%     function displayPlotPushed(~, ~)
%         % This is a callback for the Display Plot button within the display
%         % plot GUI, which will plot the information about a clicked
%         % trajectory from the drop down menu.
%
%         % Grab various fields from the TrajInfoStruct that we may need.
%         CurrentTrajectoryID = TrajInfoStruct.CurrentlySelectedTrajectory;
%         PlotDropdownMenu = TrajInfoStruct.PlotDropdownMenu;
%         TRLocal = TrajInfoStruct.TR;
%
%         % Create a new figure and plot the trajectory information as
%         % selected by the drop down menu in the display plot GUI.
%         DropdownString = PlotDropdownMenu.String;
%         DropdownValue = PlotDropdownMenu.Value;
%         DesiredPlotString = DropdownString{DropdownValue};
%         PhotonFigureWindow = findobj('Tag', 'PhotonFigure');
%         if isempty(PhotonFigureWindow)
%             PhotonFigureWindow = figure('Tag', 'PhotonFigure');
%         end
%         figure(PhotonFigureWindow); % ensure we use the correct figure
%         hold('on');
%         switch DesiredPlotString
%             case 'Photons'
%                 % Extract the photons array for this trajectory from the TR
%                 % structure, converting units if necessary.
%                 CurrentTrajBool = (cell2mat({TRLocal.TrajectoryID}) ...
%                     == CurrentTrajectoryID);
%                 FrameNum = TRLocal(CurrentTrajBool).FrameNum;
%                 FrameNum = (FrameNum-1) ...
%                     *TrajInfoStruct.UnitFlag/TrajInfoStruct.FrameRate ...
%                     + FrameNum*~TrajInfoStruct.UnitFlag;
%                 PhotonsArray = TRLocal(CurrentTrajBool).Photons;
%
%                 % Plot the Photons and label the axes.
%                 plot(FrameNum, PhotonsArray, 'x')
%                 title('Photons')
%                 xlabel(sprintf('%s (%s)', ...
%                     TrajInfoStruct.TimeDimensionString, ...
%                     TrajInfoStruct.TimeUnitString), ...
%                     'Interpreter', 'Latex')
%                 ylabel('Photons', 'Interpreter', 'Latex')
%         end
%     end
%
%     function saveMovieButtonClicked(~, ~)
%         % This is a callback function to respond to clicks of the save
%         % movie button.  Clicking that button will open a prompt to save a
%         % file so that the user can enter the filename and directory of
%         % choice.
%
%         % Obtain the filepath chosen by the user, exiting cleanly if no
%         % file was chosen.
%         [File, Path] = uiputfile({'*.mp4', 'MPEG-4 (*.mp4)'; ...
%             '*.avi', 'Uncompressed AVI (*.avi)'});
%         if (File == 0) % no file chosen
%             return
%         end
%
%         % Determine how we will be saving the video and proceed from there.
%         % NOTE: I'm basing the switch/case directly on the file extension
%         %       instead of the index available from uiputfile above because
%         %       this will still work even if I change the order of the
%         %       extensions in uiputfile.
%         [~, ~, FileExtension] = fileparts(File);
%         switch FileExtension
%             case '.mp4'
%                 % Create a VideoWriter object based on chosen filepath.
%                 VideoWriterObjectLocal = VideoWriter(...
%                     fullfile(Path, File), 'MPEG-4');
%
%                 % Set video properties specific to this format.
%                 VideoWriterObjectLocal.Quality = 100;
%             case '.avi'
%                 % Create a VideoWriter object based on chosen filepath.
%                 VideoWriterObjectLocal = VideoWriter(...
%                     fullfile(Path, File), 'Uncompressed AVI');
%         end
%         VideoWriterObjectLocal.FrameRate = TrajInfoStruct.FrameRate;
%         playMovie(VideoWriterObjectLocal);
%     end


end