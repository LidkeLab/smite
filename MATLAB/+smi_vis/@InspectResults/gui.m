function gui(obj)
%gui prepares an SR inspection GUI for the InspectResults class.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Create a figure handle for the GUI.
if (isempty(obj.GUIFigure) || ~isvalid(obj.GUIFigure))
    DefaultFigurePosition = get(0, 'defaultFigurePosition');
    obj.GUIFigure = figure('MenuBar', 'none', 'NumberTitle', 'off', ...
        'Name', 'Inspect Results', 'Tag', 'MovieGUI', ...
        'Units', 'pixels', ...
        'Position', [DefaultFigurePosition(1:2), 1000, 500]);
end

% Prepare a panel for data controls.
ControlPanelPos = [0, 0, 0.25, 1];
ControlPanel = uipanel(obj.GUIFigure, ...
    'Title', 'Control Panel', ...
    'FontWeight', 'bold', 'Units', 'normalized', ...
    'Position', ControlPanelPos);
ButtonSize = [0, 0, 0.5, 0.1];
LoadImagePos = [0, ControlPanelPos(4)-ButtonSize(4), ButtonSize(3:4)];
uicontrol('Parent', ControlPanel, ...
    'Style', 'pushbutton', 'String', 'Load Image', ...
    'Units', 'normalized', 'Position', LoadImagePos, ...
    'Callback', @loadImage);
LoadSMDPos = LoadImagePos + [ButtonSize(3), 0, 0, 0];
uicontrol('Parent', ControlPanel, ...
    'Style', 'pushbutton', 'String', 'Load SMD', ...
    'Units', 'normalized', 'Position', LoadSMDPos, ...
    'Callback', @loadSMD);
ExportButtonPos = LoadSMDPos + [0, -ButtonSize(4), 0, 0];
uicontrol('Parent', ControlPanel, ...
    'Style', 'pushbutton', 'String', 'Export SMD', ...
    'Units', 'normalized', 'Position', ExportButtonPos, ...
    'Callback', @exportSMD);
RefreshGUIPos = ButtonSize + [0, 0, 0, 0];
uicontrol('Parent', ControlPanel, ...
    'Style', 'pushbutton', 'String', 'Refresh GUI', ...
    'Units', 'normalized', 'Position', RefreshGUIPos, ...
    'Callback', @refreshGUI);

% Prepare axes for the image.
ImagePanelPos = [sum(ControlPanelPos([1, 3])), 0, 0.5, 1];
ImagePanel = uipanel(obj.GUIFigure, ...
    'Units', 'normalized', 'Position', ImagePanelPos);
obj.ImageAxes = axes(ImagePanel);
obj.ImageAxes.PositionConstraint = 'innerposition';
makeToolbar(obj.ImageAxes)

% Add a panel to contain information about selected ROI.
InfoPanelPos = [sum(ImagePanelPos([1, 3])), 0, 1-ImagePanelPos(3), 1];
InfoPanel = uipanel(obj.GUIFigure, ...
    'Title', 'ROI Information', ...
    'FontWeight', 'bold', 'Units', 'normalized', ...
    'Position', InfoPanelPos);
TextBoxSize = [0, 0, 1, 0.05];
ROITextPos = [0, InfoPanelPos(4)-TextBoxSize(4), TextBoxSize(3:4)];
ROIText = uicontrol('Parent', InfoPanel, 'Style', 'text', ...
    'String', 'ROI: []', ...
    'Units', 'normalized', ...
    'Position', ROITextPos, ...
    'HorizontalAlignment', 'left');
ROIDescriptionPos = ROITextPos + [0, -TextBoxSize(4), 0, 0];
uicontrol('Parent', InfoPanel, 'Style', 'text', ...
    'String', '        [YStart, XStart, YEnd, XEnd]', ...
    'Units', 'normalized', ...
    'Position', ROIDescriptionPos, ...
    'HorizontalAlignment', 'left');
NLocTextPos = ROIDescriptionPos + [0, -TextBoxSize(4), 0, 0];
NLocText = uicontrol('Parent', InfoPanel, 'Style', 'text', ...
    'String', '0 localizations selected', ...
    'Units', 'normalized', ...
    'Position', NLocTextPos, ...
    'HorizontalAlignment', 'left');
DatasetTextPos = NLocTextPos + [0, -TextBoxSize(4), 0, 0];
DatasetRangeText = uicontrol('Parent', InfoPanel, 'Style', 'text', ...
    'String', 'Dataset range:', ...
    'Units', 'normalized', ...
    'Position', DatasetTextPos, ...
    'HorizontalAlignment', 'left');
FrameTextPos = DatasetTextPos + [0, -TextBoxSize(4), 0, 0];
FrameRangeText = uicontrol('Parent', InfoPanel, 'Style', 'text', ...
    'String', 'Frame range:', ...
    'Units', 'normalized', ...
    'Position', FrameTextPos, ...
    'HorizontalAlignment', 'left');
DensityTextPos = FrameTextPos + [0, -TextBoxSize(4), 0, 0];
DensityText = uicontrol('Parent', InfoPanel, 'Style', 'text', ...
    'String', 'Density:', ...
    'Tooltip', 'Units of SMD pixels, not SR pixels!', ...
    'Units', 'normalized', ...
    'Position', DensityTextPos, ...
    'HorizontalAlignment', 'left');
PlotButtonPos = ButtonSize + [0, DensityTextPos(2)-ButtonSize(4), 0, 0];
uicontrol('Parent', InfoPanel, 'Style', 'pushbutton', ...
    'String', 'Display Plots', ...
    'Units', 'normalized', 'Position', PlotButtonPos, ...
    'Callback', @displayPlots);


    function updateROIInfo()
        %updateROIInfo updates the GUI to reflect changes to the ROI.
        ROIString = sprintf('%.1f, ', obj.ROI);
        ROIText.String = sprintf('ROI: [%s]', ROIString(1:(end-2)));
        if (isempty(obj.SMDIsolated) || isempty(obj.ROI))
            return
        end
        NLoc = numel(obj.SMDIsolated.FrameNum);
        NLocText.String = sprintf('%i localizations selected', NLoc);
        DatasetRangeText.String = sprintf('Dataset range: %i-%i', ...
            min(obj.SMDIsolated.DatasetNum), max(obj.SMDIsolated.DatasetNum));
        FrameRangeText.String = sprintf('Frame range: %i-%i', ...
            min(obj.SMDIsolated.FrameNum), max(obj.SMDIsolated.FrameNum));
        DensityText.String = ...
            sprintf('Density: %.4f localizations/pixel^2', ...
            NLoc / prod(obj.ROI(3:4)-obj.ROI(1:2)));
        updateAxes(obj.ImageAxes)
    end

    function displayPlots(~, ~)
        % This is a callback function to respond to clicks of the Display
        % Plots button.  Clicking that button will open a small GUI with a
        % dropdown menu to allow for display of various plots related to
        % the currently selected localizations.
        
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
        % plot GUI, which will plot the information about selected
        % localizations.
        
        % Create a new figure and plot the trajectory information as
        % selected by the drop down menu in the display plot GUI.
        DisplayPlotsFigure = findobj('Tag', 'DisplayPlotsFigure');
        if isempty(DisplayPlotsFigure)
            DisplayPlotsFigure = figure('Tag', 'DisplayPlotsFigure');
        end
        figure(DisplayPlotsFigure)
        PlotAxes = axes(DisplayPlotsFigure);
        hold(PlotAxes, 'on');
        
        % Plot the requested field from the SMD, plotting each ConnectID as
        % a different color (this might be nice for, e.g., pre-frame
        % connected data).
        DropdownString = DropdownMenu.String;
        DropdownValue = DropdownMenu.Value;
        DesiredPlotString = DropdownString{DropdownValue};
        DesiredField = obj.SMDIsolated.(DesiredPlotString);
        FrameNum = obj.SMDIsolated.FrameNum;
        DatasetNum = obj.SMDIsolated.DatasetNum;
        assert(numel(FrameNum) == numel(DesiredField), ...
            ['The selected field is not the same length as ', ...
            'SMD.FrameNum and cannot be plotted'])
        UniqueIDs = unique(obj.SMDIsolated.ConnectID);
        IDColors = lines(numel(UniqueIDs));
        for ii = 1:numel(UniqueIDs)
            CurrentBool = (obj.SMDIsolated.ConnectID == UniqueIDs(ii));
            plot(PlotAxes, FrameNum(CurrentBool) ...
                + (DatasetNum(CurrentBool)-1)*obj.SMDIsolated.NFrames, ...
                DesiredField(CurrentBool), 'x', 'Color', IDColors(ii, :))
        end
        xlabel(PlotAxes, 'Frame + (Dataset-1)*NFrames')
        ylabel(PlotAxes, DesiredPlotString)
    end

    function refreshGUI(~, ~)
        % Update the GUI to display obj.SRImage as well as any other useful
        % updates.
        cla(obj.ImageAxes)
        ImSize = size(obj.SRImage);
        imshow(obj.SRImage, [0, 1], 'Parent', obj.ImageAxes, ...
            'XData', [0.5, ImSize(2)-0.5], 'YData', [0.5, ImSize(1)-0.5])
        updateAxes(obj.ImageAxes)
        updateROIInfo()
    end

    function loadImage(Source, ~)
        %loadImage loads an SR image.

        % Load and display an image.
        try
            % Ask the user to select an image.
            Source.Enable = 'off';
            [File, Path] = uigetfile({'*.png'; '*.jpeg'}, ...
                'Select an SR reconstruction image');

            % Load the image into obj.SRImage and rescale.
            if ~isequal(File, 0)
                obj.SRImage = double(imread(fullfile(Path, File)));
                obj.SRImage = smi_vis.contrastStretch(obj.SRImage, [0, 1]);
            else
                Source.Enable = 'on';
                return
            end

            % Update the GUI.
            refreshGUI()
            Source.Enable = 'on';
        catch MException
            Source.Enable = 'on';
            rethrow(MException)
        end
    end

    function loadSMD(Source, ~)
        %loadSMD loads an SMD from a .mat file.

        % Load the SMD.
        try
            % Ask the user to select an results file.
            Source.Enable = 'off';
            [File, Path] = uigetfile('*.mat', 'Select a Results file');

            % Load the data into obj.SMD.
            if ~isequal(File, 0)
                load(fullfile(Path, File), 'SMD', 'SMF')
                obj.SMD = SMD;
                if exist('SMF', 'var')
                    % Added for back compatability with older results.
                    obj.SMF = SMF;
                end
            else
                Source.Enable = 'on';
                return
            end
        catch MException
            Source.Enable = 'on';
            rethrow(MException)
        end
        Source.Enable = 'on';
    end

    function exportSMD(~, ~)
        %exportSMD allows the user to export obj.SMDIsolated to a .mat file
        [File, Path] = uiputfile('*.mat', 'Save current ROI to a file.');
        if ~isequal(File, 0)
            SMD = obj.SMDIsolated;
            SMF = obj.SMF;
            save(fullfile(Path, File), 'SMD', 'SMF')
        end
    end

    function getROI(~, ~)
        %getROI allows selection of a rectangular ROI in the image axes.

        % Use the drawrectangle() prompt to allow user selection of the
        % ROI.
        delete(obj.ROIHandle)
        ImSize = size(obj.SRImage);
        obj.ROIHandle = drawrectangle(obj.ImageAxes, 'DrawingArea', ...
            [0, 0, ImSize(2), ImSize(1)]);

        % Define the ROI in terms of the data size (if obj.SMD is present)
        % or the image size.
        ROI = [obj.ROIHandle.Position([2, 1]), ...
            obj.ROIHandle.Position([2, 1])+obj.ROIHandle.Position([4, 3])];
        if isempty(obj.SMD)
            % In this case, only the image was provided, so we'll define
            % the ROI in image coordinates.
            obj.ROI = ROI;
        else
            ROI = ROI ./ [ImSize(1), ImSize(2), ImSize(1), ImSize(2)];
            obj.ROI = ROI .* [obj.SMD.YSize, obj.SMD.XSize, ...
                obj.SMD.YSize, obj.SMD.XSize];
        end

        % Update the ROI information panel.
        updateROIInfo()
    end

    function makeToolbar(ImageAxes)
        %makeToolbar makes a custom toolbar for the given axes.

        % Prepare a default toolbar.
        Toolbar = axtoolbar(ImageAxes, 'default');
        
        % Add a custom localization selection button.
        Icon = ones(23, 23);
        Icon(11:13, :) = 0;
        Icon(:, 11:13) = 0;
        Icon = repmat(Icon, [1, 1, 3]);
        axtoolbarbtn(Toolbar, 'push', 'Icon', Icon, ...
            'Tooltip', 'Select ROI', 'ButtonPushedFcn', @getROI);
    end

    function updateAxes(ImageAxes)
        %updateAxes updates the axes in which the SR image is displayed.

        % Ensure the toolbar is present/update it.
        makeToolbar(ImageAxes)

        % Update the axes tick labels based on obj.SMD, and if that's not
        % present, based on the size of the SR image.
        xtickformat(obj.ImageAxes, '%g')
        ytickformat(obj.ImageAxes, '%g')
        ImSize = size(obj.SRImage);
        NTicksX = numel(obj.ImageAxes.XTick);
        NTicksY = numel(obj.ImageAxes.YTick);
        obj.ImageAxes.XTick = linspace(0, ImSize(2), NTicksX);
        obj.ImageAxes.YTick = linspace(0, ImSize(1), NTicksY);
        if ~(isempty(obj.SMD) ...
                || isempty(obj.SMD.XSize) || isempty(obj.SMD.YSize))
            % In this case, we want the axes to reflect the coordinates of
            % SMD.
            XTicks = linspace(0, obj.SMD.XSize+0.5, NTicksX).';
            obj.ImageAxes.XTickLabel = num2str(XTicks, '%g');
            YTicks = linspace(0, obj.SMD.YSize+0.5, NTicksY).';
            obj.ImageAxes.YTickLabel = num2str(YTicks, '%g');
        else
            % If no SMD is present, we'll just use the image coordinates.
            obj.ImageAxes.XTickLabel = num2str(obj.ImageAxes.XTick.', '%g');
            obj.ImageAxes.YTickLabel = num2str(obj.ImageAxes.YTick.', '%g');
        end
    end


end