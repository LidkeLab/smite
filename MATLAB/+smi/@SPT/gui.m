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

% Prevent closing after a 'close' or 'close all'
GUIFigure.HandleVisibility='off';

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
ControlHandles.TestTrackButton = uicontrol(ControlPanel, ...
    'Style', 'pushbutton', 'String', 'Test Track', ...
    'FontUnits', 'normalized', 'FontSize', 0.4, ...
    'Units', 'normalized', ...
    'Position', ButtonSize + [0, ButtonSize(4), 0, 0], ...
    'callback', @testTrack);
ControlHandles.TestFitButton = uicontrol(ControlPanel, ...
    'Style', 'pushbutton', 'String', 'Test Fit', ...
    'FontUnits', 'normalized', 'FontSize', 0.4, ...
    'Units', 'normalized', ...
    'Position', ButtonSize + [0, 3*ButtonSize(4), 0, 0], ...
    'callback', @testFit);
uicontrol(ControlPanel, 'Style', 'text', 'String', 'Dataset Number:', ...
    'FontUnits', 'normalized', 'FontSize', 0.4, ...
    'Units', 'normalized', ...
    'Position', TextSize + [ButtonSize(3), 3*ButtonSize(4), 0, 0]);
ControlHandles.DatasetNumEdit = uicontrol(ControlPanel, ...
    'Style', 'edit', 'String', '1', ...
    'FontUnits', 'normalized', 'FontSize', 0.4, ...
    'Units', 'normalized', ...
    'Position', EditSize + [ButtonSize(3)+TextSize(3), 3*ButtonSize(4), 0, 0]);
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

        % Restore the SMF GUI in the panel.
        % NOTE: This is a workaround for a bug which causes detachment of
        %       the SMF GUI from the SPT class.  I haven't found the
        %       underlying cause of this bug, but this workaround seems to
        %       work okay. DJS 22/05/23
        obj.SMF = smi_core.SingleMoleculeFitting.reloadSMF(obj.SMF);
        obj.SMF.gui(SMFPanel);
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

        % Restore the SMF GUI in the panel.
        % NOTE: This is a workaround for a bug which causes detachment of
        %       the SMF GUI from the SPT class.  I haven't found the
        %       underlying cause of this bug, but this workaround seems to
        %       work okay. DJS 22/05/23
        obj.SMF = smi_core.SingleMoleculeFitting.reloadSMF(obj.SMF);
        obj.SMF.gui(SMFPanel);
    end

    function testFit(~, ~)
        % Run the smi.SMLM test fit and update some relevant class
        % properties from the results (e.g., obj.SMD).
        obj.SMLM = smi.SMLM(obj.SMF, false);
        obj.SMLM.CalledByGUI = true;
        obj.SMLM.testFit(str2double(ControlHandles.DatasetNumEdit.String));
        obj.SMD = obj.SMLM.SMD;
        obj.SMDPreThresh = obj.SMLM.SMDPreThresh;
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
