function [] = performFullAnalysis(obj)
%performFullAnalysis is the main run method for the smi.Publish class.
% This method is the main run method for the smi.Publish class, meaning
% that it can be used to perform the standard analysis expected for use of
% this class.
%
% REQUIRES:
%   matlab-instrument-control when obj.GenerateOverlayStats is true
%       (to use MIC_Reg3DTrans.findStackOffset())


% Define the results directory, which will be in the top level directory
% obj.CoverslipDir for easy access.
if isempty(obj.SaveBaseDir)
    obj.SaveBaseDir = fullfile(obj.CoverslipDir, 'Results');
end

% Define the path to the log file (a file containing info. about
% errors/warnings that happened during analysis).
StartTime = smi_helpers.genTimeString();
if isempty(obj.LogFilePath)
    FileSuffix = smi_helpers.stringMUX(...
        {obj.SMF.Data.AnalysisID, StartTime}, ...
        isempty(obj.SMF.Data.AnalysisID));
    obj.LogFilePath = fullfile(obj.SaveBaseDir, ...
        ['Log_', FileSuffix, '.mat']);
end
LogFileDir = fileparts(obj.LogFilePath);
if ~isfolder(LogFileDir)
    % The directory containing the log file doesn't exist yet, so we'll
    % make it now.
    mkdir(LogFileDir)
end
save(obj.LogFilePath, 'StartTime')

% Determine the names of the sub-directories of interest within
% obj.CoverslipDir.  These correspond to individual cells imaged during the
% experiment.
CellNamesFound = smi_helpers.getDirectoryNames(obj.CoverslipDir, 'Cell*');
CellNames = CellNamesFound;
NCellsFound = numel(CellNamesFound);
if ~isempty(obj.CellList)
    KeepCells = ones(NCellsFound, 1, 'logical');
    for ii = 1:NCellsFound
        CellNumber = sscanf(CellNames{ii}, 'Cell_%d');
        if ~ismember(CellNumber, obj.CellList)
            KeepCells(ii) = false;
        end
    end
    CellNames = CellNames(KeepCells);
end
NCells = numel(CellNames);
if (obj.Verbose > 1)
    fprintf(['Publish.performFullAnalysis(): %i ', ...
        'cell directories found:\n'], NCellsFound)
    for ii = 1:NCellsFound
        fprintf('\t%s\n', CellNamesFound{ii})
    end
    fprintf(['Publish.performFullAnalysis(): %i ', ...
        'cell directories to be analyzed:\n'], NCells)
    for ii = 1:NCells
        fprintf('\t%s\n', CellNames{ii})
    end
elseif (obj.Verbose == 1)
    fprintf(['Publish.performFullAnalysis(): analyzing %i ', ...
        'cell directories...\n'], NCells)
end

% Loop through the cell directories and analyze the contents.
for ii = 1:NCells
    % Determine if this cell should be analyzed.
    if ~isempty(obj.CellList)
        CellNumber = sscanf(CellNames{ii}, 'Cell_%d');
        if ~ismember(CellNumber, obj.CellList)
            continue
        end
    end
    if (obj.Verbose > 1)
        fprintf(['Publish.performFullAnalysis(): analyzing ', ...
            'cell directory %i of %i: %s...\n'], ii, NCells, CellNames{ii})
    end
    obj.processCell(CellNames{ii});
end

% Generate misc. stats. comparing all cells that were analyzed.
if obj.GenerateOverlayStats
    obj.genOverlayResults()
end

% Save the results.
obj.saveResults()

% Indicate completion of the analysis/generation of results.
EndTime = smi_helpers.genTimeString();
ErrorLog = obj.ErrorLog;
save(obj.LogFilePath, 'EndTime', 'ErrorLog', '-append')
if obj.Verbose
    fprintf('Results have been published to %s\n', obj.CoverslipDir);
end


end