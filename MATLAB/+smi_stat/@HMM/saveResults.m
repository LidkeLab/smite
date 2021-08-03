function saveResults(obj)
%saveResults saves useful results of the Hidden Markov model analysis.
% This method contains a collection of various plot generation codes and
% result saving codes that can be used after smi.HMM.performFullAnalysis 
% has been run.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)
    

% Make sure that the top level save directory (obj.SaveDir) exists.
if (isempty(obj.SaveDir) && ~isempty(obj.SMF.Data.ResultsDir))
    obj.SaveDir = obj.SMF.Data.ResultsDir;
else
    obj.SaveDir = pwd();
end
if ~exist(obj.SaveDir, 'dir')
    mkdir(obj.SaveDir);
end

% Extract misc. parameters from obj (to clean up the code) and define other
% useful parameters.
FrameRate = obj.SMF.Data.FrameRate;
UnitFlag = obj.UnitFlag;
NStates = numel(obj.PDFHandles);
TimeUnitString = smi_helpers.arrayMUX({'second', 'frame'}, UnitFlag);

% Create a list of directories in the top level save directory.  We will
% check if there are directories named "Condition" followed by a number,
% e.g. "Condition1", and then create the next directory following that
% pattern, e.g. "Condition2".  If the user provided a condition label,
% we won't do that.
% NOTE: We'll always just add 1 to the highest number used, so if a
%       directory named "Condition3" exists, but 1 and 2 don't, we'll still
%       just create our new directory as "Condition4".  If "Condition" is 
%       followed by a string, it shouldn't affect what we're doing here...
if isempty(obj.ConditionLabel)
    % Isolate items that follow the "Condition*" naming pattern in the
    % current directory, removing items that aren't directories just to be
    % safe.
    ConditionDirs = dir(fullfile(obj.SaveDir, 'Condition*'));
    ConditionDirNames = {ConditionDirs(...
        cell2mat({ConditionDirs.isdir})).name};
    
    % Create a suffix for our new "Condition*" directory.
    if isempty(ConditionDirNames)
        obj.ConditionLabel = 'Condition1';
    else
        % Find the maximum number is in the existing "Condition*"
        % directories.
        MaxNumber = 0; % max number found in the directory names
        for ii = 1:numel(ConditionDirNames)
            CurrentNumberString = erase(ConditionDirNames(ii), ...
                'Condition');
            MaxNumber = max(MaxNumber, str2double(CurrentNumberString{1}));
        end
        obj.ConditionLabel = sprintf('Condition%s', num2str(MaxNumber+1));
    end
end
ResultsDir = fullfile(obj.SaveDir, obj.ConditionLabel);
if ~exist(ResultsDir, 'dir')
    mkdir(ResultsDir);
end

% Isolate all pairs from the TRArray that were found to have dimerized.
TRArrayTemp = obj.TRArray;
KeepBool = zeros(size(TRArrayTemp, 2), 1);
for ii = 1:size(TRArrayTemp, 1)
    % Search for pairs that had a 1 in their StateSequence (state 1 is the
    % dimer state).
    StateSequenceCurrent = TRArrayTemp(ii, 1).StateSequence;
    KeepBool(ii) = any(StateSequenceCurrent == 1); 
end
TRArrayDimer = TRArrayTemp(logical(KeepBool), :);
TRArrayNotDimer = TRArrayTemp(~logical(KeepBool), :);
DimerDurations = cell2mat({TRArrayDimer(:, 1).DimerDurations}.');
DimerDurations = DimerDurations * (UnitFlag/FrameRate + ~UnitFlag);

% Save a .txt summary file containing the dimer->* off rates for easy
% access.
if (NStates == 3)
    KOff = obj.RateParameters(3) + obj.RateParameters(5);
    KOffSE = sqrt(obj.RateParametersSE(3)^2 + obj.RateParametersSE(5)^2);
    RateParameterKey = {'RateParameters(1): Domain -> Dimer'; ...
        'RateParameters(2): Free -> Dimer'; ...
        'RateParameters(3): Dimer -> Domain';
        'RateParameters(4): Free -> Domain';
        'RateParameters(5): Dimer -> Free';
        'RateParameters(6): Domain -> Free';
        sprintf('Rate parameter units: 1/%s', TimeUnitString)};
else
    KOff = obj.RateParameters(2);
    KOffSE = obj.RateParametersSE(2);
    RateParameterKey = {'RateParameters(1): Free -> Dimer'; ...
        'RateParameters(2): Dimer -> Free'; ...
        sprintf('Rate parameter units: 1/%s', TimeUnitString)};
end
KOff = KOff * (UnitFlag*FrameRate + ~UnitFlag); % convert to desired units
KOffSE = KOffSE * (UnitFlag*FrameRate + ~UnitFlag);
FileID = fopen(fullfile(obj.SaveDir, 'OffRateSummary.txt'), 'a');
fprintf(FileID, '%s: k_off = (%0.4f+-%0.4f) / %s, %i dimer events\n', ...
    obj.ConditionLabel, KOff, KOffSE, TimeUnitString, ...
    numel(DimerDurations));
fclose(FileID);

% Save the rate parameters in the top level directory, appending the
% filename and variable names with the specified condition label.
RateParameterNameString = sprintf('RateParameters%s', ...
    obj.ConditionLabel);
RateParametersTemp = obj.RateParameters ...
    * (UnitFlag*FrameRate + ~UnitFlag); % convert to desired units
eval(sprintf('%s = RateParametersTemp;', ...
    RateParameterNameString)); % create local var. with descriptive name
RateParameterSENameString = sprintf('RateParametersSE%s', ...
    obj.ConditionLabel);
RateParametersSETemp = obj.RateParametersSE ...
    * (UnitFlag*FrameRate + ~UnitFlag); % convert to desired units
eval(sprintf('%s = RateParametersSETemp;', ...
    RateParameterSENameString)); % create local var. with descriptive name
save(fullfile(obj.SaveDir, [RateParameterNameString, '.mat']), ...
    RateParameterNameString, RateParameterSENameString, ...
    'RateParameterKey');

% Save a copy of this class instance in the ResultsDir.
HMMNameString = 'HMM';
eval(sprintf('%s = obj;', HMMNameString)); % create local variable
save(fullfile(ResultsDir, ...
    sprintf('%s%s.mat', HMMNameString, obj.ConditionLabel)), ...
    HMMNameString, '-v7.3');

% Save these sub-TRArrays in the appropriate directories (for easy user
% access).
DimerDir = fullfile(ResultsDir, 'Dimers');
if ~exist(DimerDir, 'dir')
    mkdir(DimerDir);
end
NotDimerDir = fullfile(ResultsDir, 'Other');
if ~exist(NotDimerDir, 'dir')
    mkdir(NotDimerDir);
end
save(fullfile(DimerDir, 'TRArrayDimer.mat'), 'TRArrayDimer', '-v7.3');
save(fullfile(NotDimerDir, 'TRArrayNotDimer.mat'), 'TRArrayNotDimer', ...
    '-v7.3');

% Generate some interesting plots (if desired).
if obj.GeneratePlots
    % Create a histogram of the event durations found from the Viterbi
    % algorithm, with the exponential distribution found in HMM plotted on
    % top.
    FigureHandle = figure('Visible', 'off');
    PlotAxes = axes(FigureHandle);
    hold(PlotAxes, 'on');
    BinWidth = max(1, round(mean(DimerDurations)));
    [~, BinEdges] = histcounts(DimerDurations, 'BinWidth', BinWidth);
    histogram(PlotAxes, DimerDurations, 'BinWidth', BinWidth, ...
        'Normalization', 'pdf')
    DurationRange = linspace(BinEdges(1), BinEdges(end), 1e3);
    ExpFit = KOff * exp(-KOff*DurationRange);
    plot(PlotAxes, DurationRange, ExpFit);
    legend(PlotAxes, {'Observed durations', ...
        sprintf('k_{off} = %0.4f / %s', KOff, TimeUnitString)}, ...
        'Location', 'best');
    title(PlotAxes, sprintf('%s, %i dimer events', ...
        obj.ConditionLabel, numel(DimerDurations)), 'Interpreter', 'none')
    xlabel(PlotAxes, sprintf('Dimer duration (%ss)', TimeUnitString), ...
            'Interpreter', 'Latex')
    ylabel(PlotAxes, 'PDF', 'Interpreter', 'Latex')
    FileName = fullfile(DimerDir, 'DimerDurationsHistogram');
    saveas(FigureHandle, FileName, 'png');
    close(FigureHandle);
    
    % Create a single plot with several interesting subplots that can be
    % viewed quickly as a summary.
    FigureHandle = figure('Units', 'inches', ...
        'Position', [0, 0, 8.5, 11], ...
        'PaperSize', [8.5, 11]);
    DisplayParams.ChannelNames = obj.ChannelNames;
    DisplayParams.StateNames = obj.StateNames;
    DisplayParams.StateColormap = colormap(lines(NStates));
    DisplayParams.UnitFlag = UnitFlag;
    DisplayParams.MinXYRange = 10; % pixels
    DisplayParams.MaxYDisplaySep = obj.MaxSeparation; % pixels
    for ii = 1:size(TRArrayDimer, 2)
        % Generate the plot.
        DisplayParams.PairNumber = ii;
        obj.createSummaryPlot(TRArrayDimer(:, ii), DisplayParams, ...
            FigureHandle);
        
        % Save the plot.
        FileName = fullfile(DimerDir, sprintf('DimerPair%i', ii));
        saveas(FigureHandle, FileName, 'fig');
        saveas(FigureHandle, FileName, 'png');
        clf(FigureHandle);
    end
    close(FigureHandle);
end

% Generate dimer movies (if requested).
if obj.GenerateMovies
    DisplayParams.ChannelNames = obj.ChannelNames;
    DisplayParams.UnitFlag = 1;
    DisplayParams.MaxTrajDisplayLength = inf;
    DisplayParams.AutoClip = 1;
    DisplayParams.MinXYRange = 20; % pixels
    DisplayParams.NPadPixels = 5; % pixels
    DisplayParams.NPadFrames = 10; % frames
    DisplayParams.IndicateDimer = 1;
    DisplayParams.IndicateHMMData = 1;
    obj.createAllMovies(TRArrayDimer, DimerDir, '', DisplayParams);
end


end