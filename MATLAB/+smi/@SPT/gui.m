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
ControlHandles.Track = uicontrol(ControlPanel, ...
    'Style', 'pushbutton', 'String', 'Track', ...
    'FontUnits', 'normalized', 'FontSize', 0.4, ...
    'Units', 'normalized', ...
    'Position', ButtonSize, ...
    'callback', @track);
ControlHandles.TestFitButton = uicontrol(ControlPanel, ...
    'Style', 'pushbutton', 'String', 'Test Track', ...
    'FontUnits', 'normalized', 'FontSize', 0.4, ...
    'Units', 'normalized', ...
    'Position', ButtonSize + [0, ButtonSize(4), 0, 0], ...
    'callback', @testTrack);
ControlHandles.MakeMovie = uicontrol(ControlPanel, ...
    'Style', 'pushbutton', 'String', 'Movie GUI', ...
    'FontUnits', 'normalized', 'FontSize', 0.4, ...
    'Units', 'normalized', ...
    'Position', ButtonSize + [ControlPanel.Position(3)-ButtonSize(3), 0, 0, 0], ...
    'callback', @makeMovieGUI);

    function track(Source, ~)
        % Track the data based on the SMF GUI parameters.
        try
            Source.Enable = 'off';
            FindFilesInit = obj.FindFiles;
            obj.FindFiles = false;
            if (numel(obj.SMF.Data.FileName) > 1)
                % If multiple filenames are present, we'll try
                % batch-tracking.
                obj.batchTrack();
            else
                obj.performFullAnalysis()
            end
            Source.Enable = 'on';
            obj.FindFiles = FindFilesInit;
        catch MException
            Source.Enable = 'on';
            obj.FindFiles = FindFilesInit;
            rethrow(MException)
        end
    end

    function testTrack(Source, ~)
        % Track the data based on the SMF GUI parameters.
        TestFlagInit = obj.IsTestRun;
        VerboseInit = obj.Verbose;
        FindFilesInit = obj.FindFiles;
        obj.IsTestRun = true;
        obj.Verbose = 3;
        obj.FindFiles = false;
        try
            Source.Enable = 'off';
            if (numel(obj.SMF.Data.FileName) > 1)
                % For multiple files, we'll dispatch on batchTrack().
                obj.batchTrack();
            else
                obj.performFullAnalysis()
            end
            Source.Enable = 'on';
        catch MException
            Source.Enable = 'on';
            obj.IsTestRun = TestFlagInit;
            obj.Verbose = VerboseInit;
            obj.FindFiles = FindFilesInit;
            rethrow(MException)
        end
        obj.IsTestRun = TestFlagInit;
        obj.Verbose = VerboseInit;
        obj.FindFiles = FindFilesInit;
    end

    function makeMovieGUI(~, ~)
        % Prepare the movie maker.
        MovieMaker = smi_vis.GenerateMovies;
        MovieMaker.TR = obj.TR;
        MovieMaker.SMD = obj.SMD;
        MovieMaker.RawData = obj.ScaledData;
        MovieMaker.SMF = copy(obj.SMF);
        MovieMaker.gui()
    end


end