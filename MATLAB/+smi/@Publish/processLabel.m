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
        fprintf('\t\t%s\n', DataFileNames{ii})
    end
end

% Load and analyze the data for the current label, looping through datasets
% if needed.
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
    SaveDir = fullfile(obj.SaveBaseDir, CellName, LabelName);
    
    % Generate figures associated with the brightfield registration of the
    % cell (if that data exists).
    if obj.GenerateImagingStats
        % Define the path to the .h5 data file.
        FilePath = fullfile(DataDir, DataFileNames{ii});
        
        % Create a list of the groups available in each dataset.
        H5FileStruct = h5info(FilePath);
        FileGroupList = ...
            {H5FileStruct.Groups.Groups.Groups(1).Groups.Name};
        if any(contains(FileGroupList, 'AlignReg'))
            if (obj.Verbose > 1)
                fprintf(['\t\tPublish.processLabel(): ', ...
                    'Generating cell registration results...\n'])
            end
            [~, FileNameNoExt] = fileparts(DataFileNames{ii});
            ImagingStatsSaveDir = fullfile(SaveDir, FileNameNoExt);
            [AlignResultsStruct] = ...
                obj.genAlignResults(FilePath, ImagingStatsSaveDir);
            
            % Store information from the AlignResultsStruct into a field
            % associated with the current dataset in the
            % PublishedResultsStruct.
            FieldNames = fieldnames(AlignResultsStruct);
            LabelNumber = str2double(regexp(LabelName, '\d\d', 'match'));
            CellNumber = str2double(regexp(CellName, '\d\d', 'match'));
            obj.ResultsStruct(CellNumber, LabelNumber).Cell = CellName;
            obj.ResultsStruct(CellNumber, LabelNumber).Label = LabelName;
            for jj = 1:numel(FieldNames)
                obj.ResultsStruct(CellNumber, LabelNumber). ...
                    (FieldNames{jj}) = AlignResultsStruct.(FieldNames{jj});
            end
        end
    end
    
    % Generate the super-resolution results using the smi.SMLM class.
    obj.SMLM.SMF.Data.FileName = DataFileNames(ii);
    if obj.GenerateSR
        if (obj.Verbose > 1)
            fprintf(['\t\tPublish.processLabel(): ', ...
                'Analyzing file %s...\n'], ...
                fullfile(DataDir, DataFileNames{ii}))
        end
        
        % Place this in a try/catch so that we can still proceed with
        % the other analyses if this fails.
        try
            % Perform SR analysis, but don't save the results yet!
            obj.SMLM.SMF.Data.ResultsDir = SaveDir;
            obj.SMLM.FullvsTest = true;
            obj.SMLM.analyzeAll();

            % Re-shift the XY coordinates w.r.t. the best registration
            % dataset.
            if (obj.ShiftToReg && obj.GenerateImagingStats)
                % Load the brightfield images taken just before each
                % sequence at the focal plane.
                FocusImages = ...
                    smi_core.LoadData.readH5File(FilePath, 'FocusImages');
                
                % Determine which focus image has the least shift w.r.t.
                % the reference (in XY, ignoring Z).
                RefImage = AlignResultsStruct.AlignReg(1).Data.Image_Reference;
                obj.SMLM.SMD = obj.shiftToBestReg(obj.SMLM.SMD, ...
                    RefImage, FocusImages);
            end
            
            % Save the SR results.
            obj.SMLM.saveResults();
            
            % Copy the results into a more accessible directory (it's nice
            % to have them all in one place for certain analyses).
            [~, FileName] = fileparts(DataFileNames{ii});
            ResultsFile = dir(fullfile(obj.SMLM.SMF.Data.ResultsDir, ...
                FileName, '*Results.mat'));
            if ~isempty(ResultsFile)
                ResultsStructDir = ...
                    fullfile(obj.SaveBaseDir, 'ResultsStructs');
                if ~isfolder(ResultsStructDir)
                    mkdir(ResultsStructDir)
                end
                NewName = [CellName, '_', LabelName, '_', 'Results.mat'];
                copyfile(fullfile(ResultsFile.folder, ResultsFile.name), ...
                    fullfile(ResultsStructDir, NewName));
            end
        catch MException
            if obj.Verbose
                warning(['Publish.processLabel(): ', ...
                    'Unable to generate SR images for %s\n%s, %s'], ...
                    DataFileNames{ii}, ...
                    MException.identifier, MException.message)
            end
            
            % Store the error information in the log file.
            obj.ErrorLog = [obj.ErrorLog; ...
                {CellName, LabelName, MException}];
            ErrorLog = obj.ErrorLog;
            save(obj.LogFilePath, 'ErrorLog', '-append')
        end
    end
end


end