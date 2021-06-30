function gui(obj)
%gui prepares a movie generation GUI for the GenerateMovies class.
%
% INPUTS:
%   GUIParent: The 'Parent' of this GUI, e.g., a figure handle.
%              (Default = figure(...))

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Create a figure handle for the GUI.
PlotFigure = figure('MenuBar', 'none', ...
    'Name', 'Movie maker', 'NumberTitle', 'off', ...
    'Units', 'pixels');

% Define some useful parameters.
if obj.Params.UnitFlag
    % UnitFlag == 1 corresponds to physical units (micrometers and seconds)
    LengthUnitString = '$\mu m$';
    TimeDimensionString = 'Time';
    TimeUnitString = 's';
else
    % UnitFlag == 0 corresponds to camera units (pixels and frames).
    LengthUnitString = 'pixels';
    TimeDimensionString = 'Frame';
    TimeUnitString = 'frames';
end

% Prepare axes for the movie.
AxesPanelPos = [0, 0.11, 0.71, 0.89];
AxesPanel = uipanel(PlotFigure, ...
    'Units', 'normalized', 'Position', AxesPanelPos);
obj.MovieAxes = axes(AxesPanel);
obj.MovieAxes.ActivePositionProperty = 'position';

% Add some controls for the movie.
MovieControlPanel = uipanel(PlotFigure, 'Title', 'Movie Controls', ...
    'FontWeight', 'bold', 'Units', 'normalized', ...
    'Position', AxesPanelPos + [0, -0.11, 0, -0.78]);
FrameSliderPos = [0.1, 0.5, 0.89, 0.5];
FrameSlider = uicontrol('Parent', MovieControlPanel, ...
    'Style', 'slider', ...
    'Units', 'normalized', 'Position', FrameSliderPos, ...
    'HorizontalAlignment', 'left', ...
    'Min', 0, 'Max', 1, ...
    'Enable', 'off', 'Callback', @frameSlider);
PlayButtonPos = [0, FrameSliderPos(2), ...
    0.9 * (1-FrameSliderPos(3)), FrameSliderPos(4)];
PlayButton = uicontrol('Parent', MovieControlPanel, ...
    'Style', 'pushbutton', 'String', 'PLAY', ...
    'Units', 'normalized', 'Position', PlayButtonPos, ...
    'Enable', 'off', 'Callback', @playMovieCallback);
FrameNumDisplayPos = [FrameSliderPos(1), 0, ...
    FrameSliderPos(3), 1 - FrameSliderPos(4)];
FrameNumDisplay = uicontrol('Parent', MovieControlPanel, ...
    'Style', 'text', ...
    'Units', 'normalized', 'Position', FrameNumDisplayPos, ...
    'HorizontalAlignment', 'left');

% Add a panel to the figure to contain trajectory information about clicked
% trajectories.
TrajInfoPanelPos = ...
    [AxesPanelPos(1)+AxesPanelPos(3), AxesPanelPos(2)+AxesPanelPos(4)/2, ...
    1-AxesPanelPos(3), AxesPanelPos(4)/2];
TrajInfoPanel = uipanel(PlotFigure, ...
    'Title', 'Trajectory Information', ...
    'FontWeight', 'bold', 'Units', 'normalized', ...
    'Position', TrajInfoPanelPos);
TrajInfoTitlePos = [0, 0.8, 1, 0.2];
uicontrol('Parent', TrajInfoPanel, 'Style', 'text', ...
    'Units', 'normalized', 'Position', TrajInfoTitlePos, ...
    'String', sprintf(...
    '(click a trajectory to display trajectory information)'));
TrajInfoTextPos = [0, TrajInfoTitlePos(2)-TrajInfoTitlePos(4), 0.5, 0.1];
uicontrol('Parent', TrajInfoPanel, 'Style', 'text', ...
    'String', 'Trajectory ID: ', ...
    'Units', 'normalized', ...
    'Position', TrajInfoTextPos, ...
    'HorizontalAlignment', 'left');
uicontrol('Parent', TrajInfoPanel, 'Style', 'edit', ...
    'Units', 'normalized', ...
    'Position', TrajInfoTextPos ...
    + [TrajInfoTextPos(3), 0, 0, 0], ...
    'HorizontalAlignment', 'left', ...
    'Callback', @trajectorySelectedEdit, 'Tag', 'TrajIDDisplay');
uicontrol('Parent', TrajInfoPanel, 'Style', 'text', ...
    'String', sprintf('Start %s (%s):', ...
    TimeDimensionString, TimeUnitString), ...
    'Units', 'normalized', ...
    'Position', TrajInfoTextPos ...
    + [0, -TrajInfoTextPos(4), ...
    (1-TrajInfoTextPos(3)), 0], ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'StartFrameDisplay');
uicontrol('Parent', TrajInfoPanel, 'Style', 'text', ...
    'String', sprintf('End %s (%s):', ...
    TimeDimensionString, TimeUnitString), ...
    'Units', 'normalized', 'Position', TrajInfoTextPos ...
    + [0, -2*TrajInfoTextPos(4), ...
    (1-TrajInfoTextPos(3)), 0], ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'EndFrameDisplay');
uicontrol('Parent', TrajInfoPanel, 'Style', 'pushbutton', ...
    'String', 'Display Plots', ...
    'Units', 'normalized', 'Position', [0, 0, 1, 0.2], ...
    'Callback', @displayPlots);

% Add a panel to the figure to contain save information.
SaveMoviePanelPos = [TrajInfoPanelPos(1), AxesPanelPos(2), ...
    TrajInfoPanelPos(3), AxesPanelPos(4)-TrajInfoPanelPos(4)];
SaveMoviePanel = uipanel(PlotFigure, 'Title', 'Save Movie', ...
    'FontWeight', 'bold', 'Units', 'normalized', ...
    'Position', SaveMoviePanelPos);
uicontrol('Parent', SaveMoviePanel, 'Style', 'pushbutton', ...
    'String', 'Save Movie', ...
    'Units', 'normalized', 'Position', [0, 0, 1, 0.2], ...
    'Callback', @saveMovieButtonClicked);

%
%     function frameSlider(Source, ~)
%         % This is a callback function for the event that the user has slid
%         % the FrameSlider slidebar.  This method will determine the
%         % location of the slidebar, convert this to a frame within the
%         % movie, and then plot the trajectories of the movie as they exists
%         % up through that frame.
%
%         % Get the location of the slide bar, as well as the min and max
%         % values possible for the slider.
%         SliderValue = get(Source, 'Value');
%         SliderMin = get(Source, 'Min');
%         SliderMax = get(Source, 'Max');
%
%         % Convert the SliderValue to an appropriate frame within the movie.
%         % NOTE: Since the frame selected with the slidebar may not
%         %       necessarily exist in the movie, we instead select the frame
%         %       in the movie that is nearest to that chosen by the slider.
%         MinFrame = min(TrajInfoStruct.MovieFrames);
%         MaxFrame = max(TrajInfoStruct.MovieFrames);
%         SelectedFrame = round((1 / (SliderMax-SliderMin)) ...
%             * ((MaxFrame-MinFrame)*SliderValue ...
%             + SliderMax*MinFrame - SliderMin*MaxFrame)); % ideal endframe
%         [~, NearestIndex] = min(...
%             abs(SelectedFrame - TrajInfoStruct.MovieFrames));
%         EndFrame = TrajInfoStruct.MovieFrames(NearestIndex);
%
%         % Prepare the axes for the plots.
%         cla(TrajInfoStruct.PlotAxes);
%         iptsetpref('ImshowAxesVisible', 'on');
%         hold(TrajInfoStruct.PlotAxes, 'on');
%
%         % Re-plot the trajectories up through FrameNumber, plotting the
%         % entire trajectory if EndFrame is the last frame of the movie
%         % (i.e., don't implement a MaxTrajDisplayLength for the last frame)
%         if EndFrame == max(TrajInfoStruct.MovieFrames)
%             % The chosen EndFrame is the last frame of the movie.
%             plotTrajectories(min(TrajInfoStruct.MovieFrames), ...
%                 EndFrame, inf);
%         else
%             plotTrajectories(min(TrajInfoStruct.MovieFrames), EndFrame, ...
%                 TrajInfoStruct.MaxTrajDisplayLength);
%         end
%
%         % Update the frame number display within the GUI to show the
%         % currently selected frame of the movie.
%         TrajInfoStruct.FrameNumDisplay.String = sprintf(...
%             'Frame %s of %s', num2str(EndFrame), ...
%             num2str(max(TrajInfoStruct.MovieFrames)));
%     end
%
%     function playMovieCallback(Source, ~)
%         % This is a callback function for the PLAY button, which will call
%         % the function playMovie to replay the movie without needing to use
%         % the slidebar.
%
%         % Disable the PLAY button so it can't be clicked multiple times.
%         Source.Enable = 'off';
%
%         % Clear the movie axes to ensure old plots aren't still remaining.
%         cla(TrajInfoStruct.PlotAxes);
%
%         % Play the movie.
%         playMovie();
%
%         % Re-enable the PLAY button.
%         Source.Enable = 'on';
%     end
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
%     function trajectoryClicked(~, ~, ii)
%         % This is a callback function for the event that the user has
%         % clicked on a trajectory (with trajectory ID input as ii) within
%         % the movie.
%
%         % Pass the trajectory ID ii along to the trajectorySelected
%         % function.
%         trajectorySelected(ii);
%     end
%
%     function trajectorySelected(ii)
%         % This is a function used to display misc. information about a
%         % specific trajectory within the movie figure.
%
%         % Remove previous indicator(s) that a trajectory had been selected
%         % prior to this event, ensuring that the line is still valid before
%         % attempting to modify it.
%         [TrajInfoStruct.LineHandles(TrajInfoStruct.LineHandles.isvalid ...
%             & isgraphics(TrajInfoStruct.LineHandles)).LineWidth] = ...
%             deal(0.5);
%
%         % Add indicators to the selected trajectory to emphasize which one
%         % was selected.
%         TrajInfoStruct.LineHandles(...
%             TrajInfoStruct.LineHandleIDMap == ii).LineWidth = 3;
%
%         % Update the CurrentlySelectedTrajectory field within the
%         % TrajInfoStruct so that other callbacks can access this
%         % information.
%         TrajInfoStruct.CurrentlySelectedTrajectory = ii;
%
%         % Grab useful information about the trajectory to be displayed,
%         % ensuring the units are consistent with the movie.
%         CurrentTrajBool = (cell2mat({TrajInfoStruct.TR.TrajectoryID}) ...
%             == ii);
%         FrameNum = TrajInfoStruct.TR(CurrentTrajBool).FrameNum;
%         FrameNum = (FrameNum-1) ...
%             * TrajInfoStruct.UnitFlag/TrajInfoStruct.FrameRate ...
%             + FrameNum*~TrajInfoStruct.UnitFlag;
%         StartFrame = min(FrameNum); % first frame trajectory appears in
%         EndFrame = max(FrameNum); % last frame trajectory appears in
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
%                         TrajInfoStruct.TimeUnitString, StartFrame);
%                 case 'EndFrameDisplay'
%                     TrajInfoStruct.TrajInfoPanel.Children(kk).String = ...
%                         sprintf('End %s (%s): %i', ...
%                         TrajInfoStruct.TimeDimensionString, ...
%                         TrajInfoStruct.TimeUnitString, EndFrame);
%             end
%         end
%     end
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