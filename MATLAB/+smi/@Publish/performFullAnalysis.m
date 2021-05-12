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

% Determine the names of the sub-directories of interest within
% obj.CoverslipDir.  These correspond to individual cells imaged during the
% experiment.
CellNames = smi_helpers.getDirectoryNames(obj.CoverslipDir, 'Cell*');
NCells = numel(CellNames);
if (obj.Verbose > 1)
    fprintf(['Publish.performFullAnalysis(): %i ', ...
        'cell directories found:\n'], NCells)
    for ii = 1:NCells
        fprintf('\t%s\n', CellNames{ii})
    end
elseif (obj.Verbose == 1)
    fprintf(['Publish.performFullAnalysis(): analyzing %i ', ...
        'cell directories...\n'], NCells)
end

% Loop through the cell directories and analyze the contents.
for ii = 1:NCells
    if (obj.Verbose > 1)
        fprintf(['Publish.performFullAnalysis(): analyzing ', ...
            'cell %i of %i...\n'], ii, NCells)
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
if obj.Verbose
    fprintf('Results have been published to %s\n', obj.CoverslipDir);
end


end