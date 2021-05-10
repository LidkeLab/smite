function processLabel(obj, CellName, LabelName)
%processLabel processes the SR data for the specified label.
% This method will find the data in
%
% INPUTS:
%   CellName: Char. array/string of the Cell* directory (i.e., the distinct
%             Cell* names in the directory structure
%             *\Cell*\Label*\Data*.h5)

% Created by:
%   David J. Schodt (Lidke Lab 2021)


% Specify the directory containing this labels data and create a list of
% the data files.
DataDir = fullfile(obj.CoverslipDir, CellName, LabelName);
obj.SMF.Data.FileDir = DataDir;
DataFileStruct = dir(fullfile(DataDir, 'Data*'));
DataFileNames = {DataFileStruct.name};
NDataFiles = numel(DataFileNames);
if (obj.Verbose > 1)
    fprintf('\tPublish.processLabel(): %i files found:\n', NDataFiles)
    for ii = 1:NDataFiles
        fprintf('\t%s\n', DataFileNames{ii})
    end
end

% Load and analyze the data for the current label, looping through datasets
% if needed.
SMLM = smi.SMLM(obj.SMF);
for ii = 1:NDataFiles
    % Determine if this was a bleaching round and decide if we should
    % analyze it.
    if (contains(DataFileNames{ii}, 'bleach') && ~obj.AnalyzeBleaching)
        continue
    end
    
    % Display some message to the Command Window to show progression
    % through the workflow.
    if obj.Verbose
        fprintf('\tPublish.processLabel(): Analyzing %s...\n', ...
            fullfile(DataDir, DataFileNames{ii}))
    end
    
    % Define the save directory for this specific file.
    [~, FileNameNoExtension] = fileparts(DataFileNames{ii});
    SaveDir = fullfile(obj.SaveBaseDir, CellName, LabelName, ...
        FileNameNoExtension);
    
    % Generate the super-resolution results using the smi.SMLM class.
    SMLM.SMF.Data.FileName = DataFileNames(ii);
    if obj.GenerateSR
        SMLM.SMF.Data.ResultsDir = SaveDir;
        if (obj.Verbose > 1)
            fprintf(['\t\tPublish.processLabel(): ', ...
                'Analyzing file %s...\n'], ...
                fullfile(DataDir, DataFileNames{ii}))
        end
        try
            % Place this in a try/catch so that we can still proceed with
            % the other analyses if this fails.
            SMLM.analyzeAll()
        catch MException
            if obj.Verbose
                warning(['Publish.processLabel(): ', ...
                    'Unable to generate SR images for %s\n%s, %s'], ...
                    DataFileNames{ii}, ...
                    MException.identifier, MException.message)
            end
        end
    end
    
    % Generate figures associated with the brightfield registration of the
    % cell (if that data exists).
    if obj.GenerateImagingStats
        % Define the path to the .h5 data file.
        FilePath = fullfile(DataDir, DataFileNames{ii});
        
        % Create a list of the groups available in each dataset
        H5FileStruct = h5info(FilePath);
        FileGroupList = ...
            {H5FileStruct.Groups.Groups.Groups(1).Groups.Name};
        if any(contains(FileGroupList, 'AlignReg'))
            if (obj.Verbose > 1)
                fprintf(['\t\tPublish.processLabel(): ', ...
                    'Generating cell registration results...\n'])
            end
            [AlignResultsStruct] = ...
                obj.genAlignResults(FilePath, SaveDir);
            
            % Store information from the AlignResultsStruct into a field
            % associated with the current dataset in the
            % PublishedResultsStruct.
            LabelField = sprintf('%s_%s', CellName, LabelName);
            DataField = strrep(DataFileNameNoExtension, '-', '_');
            FieldNames = fieldnames(AlignResultsStruct);
            for jj = 1:numel(FieldNames)
                obj.PublishedResultsStruct.(LabelField). ...
                    (DataField).(FieldNames{jj}) = ...
                    AlignResultsStruct.(FieldNames{jj});
            end
        end
    end
end


end