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

% Grab some useful information from the file name.
CellNumber = regexp(CellName, '\d*', 'match');
CellNumber = str2double(CellNumber{1});
LabelNumber = regexp(LabelName, '\d*', 'match');
LabelNumber = str2double(LabelNumber{1});

% Load and analyze the data for the current label, looping through datasets
% if needed.
for ii = 1:NDataFiles
    % Determine if this was a bleaching round and decide if we should
    % analyze it.
    if (contains(DataFileNames{ii}, 'bleach') && ~obj.AnalyzeBleaching)
        continue
    end

    % Load the focus images taken before/after each sequence (these may be
    % used later).
    % NOTE: This only works correctly if we have one file per label 
    %       (excluding 'bleaching' files).
    FilePath = fullfile(DataDir, DataFileNames{ii});
    H5FileStruct = h5info(FilePath);
    FileGroupList = {H5FileStruct.Groups.Groups.Groups(1).Groups.Name};
    FocusImagesPresent = any(contains(FileGroupList, 'FocusImages'));
    if FocusImagesPresent
        % If the FocusImages field is present, we'll use those.
        FocusImageStruct = smi_core.LoadData.readH5File(...
            FilePath, 'FocusImages');
        FocusImages = cell(numel(FocusImageStruct), 1);
        for nn = 1:numel(FocusImageStruct)
            FocusImages{nn} = FocusImageStruct(nn).Data.PreSeqImages;
        end
    else
        % If the FocusImages field isn't present (older data), we can still
        % try to make use of the last brightfield image taken before each
        % sequence.  The AlignReg structure has a similar organization, so
        % we can just add a new field and store it in
        % obj.FocusImageStructs.
        AlignReg = smi_core.LoadData.readH5File(FilePath, 'AlignReg');
        for nn = 1:numel(AlignReg)
            FocusImages{nn} = AlignReg(nn).Data.Image_Current;
            AlignReg(nn).Data.PreSeqImages = AlignReg(nn).Data.Image_Current;
        end
        FocusImageStruct = AlignReg;
    end
    obj.FocusImageStructs{CellNumber, LabelNumber} = FocusImageStruct;

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
    H5FileStruct = h5info(FilePath);
    FileGroupList = ...
        {H5FileStruct.Groups.Groups.Groups(1).Groups.Name};
    if obj.GenerateImagingStats
        % Generate the results.
        if any(contains(FileGroupList, 'AlignReg'))
            if (obj.Verbose > 1)
                fprintf(['\t\tPublish.processLabel(): ', ...
                    'Generating cell registration results...\n'])
            end
            [~, FileNameNoExt] = fileparts(DataFileNames{ii});
            ImagingStatsSaveDir = fullfile(SaveDir, FileNameNoExt);
            [AlignResultsStruct] = ...
                obj.genAlignResults(FilePath, ImagingStatsSaveDir);
            AlignReg = AlignResultsStruct.AlignReg; % grab for later.
            obj.AlignRegStructs{CellNumber, LabelNumber} = AlignReg;

            % Store information from the AlignResultsStruct into a field
            % associated with the current dataset in the
            % PublishedResultsStruct.
            FieldNames = fieldnames(AlignResultsStruct);
            obj.ResultsStruct(CellNumber, LabelNumber).Cell = CellName;
            obj.ResultsStruct(CellNumber, LabelNumber).Label = LabelName;
            for jj = 1:numel(FieldNames)
                obj.ResultsStruct(CellNumber, LabelNumber). ...
                    (FieldNames{jj}) = AlignResultsStruct.(FieldNames{jj});
            end
        end
    end

    % Generate the super-resolution results using the smi.SMLM class.
    obj.SMLM.SMF.Data.FileDir = DataDir;
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
            if (obj.UseBrightfieldDC && FocusImagesPresent)
                % If we're using brightfield DC, we need to apply it before
                % using the post-processing DC algorithm.
                obj.SMLM.SMF.DriftCorrection.On = false;
                obj.SMLM.analyzeAll();

                % Perform the brightfield DC.
                Params.CorrParams.SuppressWarnings = true;
                AlignReg = smi_core.LoadData.readH5File(FilePath, 'AlignReg');
                obj.AlignRegStructs{CellNumber, LabelNumber} = AlignReg;
                obj.SMLM.SMD = ...
                    smi_core.DriftCorrection.driftCorrectBF(...
                    obj.SMLM.SMD, obj.SMLM.SMF, ...
                    AlignReg(1).Data.Image_Reference, ...
                    obj.FocusImageStructs{CellNumber, LabelNumber}, Params);

                % Run the post-processing drift correction if needed.
                if obj.SMF.DriftCorrection.On
                    DC = smi_core.DriftCorrection;
                    DC.SMF = obj.SMF;
                    DC.BFRegistration = true;
                    obj.SMLM.SMD = DC.driftCorrectKNN(obj.SMLM.SMD);
                end
            else
                if ((obj.Verbose>0) && ~FocusImagesPresent)
                    warning(['Brightfield drift-correction cannot be ', ...
                        'applied: FocusImages not present in .h5 file'])
                end
                obj.SMLM.analyzeAll();
            end

            % Re-shift the XY coordinates w.r.t. the best registration
            % dataset.
            if obj.ShiftToReg
                % Check if the alignment results are available, loading
                % them if not.
                if (~exist('AlignReg', 'var') || isempty(AlignReg))
                    AlignReg = smi_core.LoadData.readH5File(FilePath, 'AlignReg');
                    obj.AlignRegStructs{CellNumber, LabelNumber} = AlignReg;
                end

                % Determine which focus image has the least shift w.r.t.
                % the reference (in XY, ignoring Z).
                RefImage = AlignReg(1).Data.Image_Reference;
                obj.SMLM.SMD = obj.shiftToBestReg(obj.SMLM.SMD, ...
                    RefImage, FocusImages);
            end

            % Save the SR results.
            obj.SMLM.saveResults();

            % Copy the results into a more accessible directory (it's nice
            % to have them all in one place for certain analyses).
            [~, FileName] = fileparts(DataFileNames{ii});
            FileNameSaved = [FileName, ...
                smi_helpers.arrayMUX({'_', ''}, ...
                isempty(obj.SMLM.SMF.Data.AnalysisID)), ...
                obj.SMLM.SMF.Data.AnalysisID];
            ResultsFile = dir(fullfile(obj.SMLM.SMF.Data.ResultsDir, ...
                FileNameSaved, '*Results.mat'));
            if ~isempty(ResultsFile)
                ResultsStructDir = ...
                    fullfile(obj.SaveBaseDir, 'ResultsStructs');
                if ~isfolder(ResultsStructDir)
                    mkdir(ResultsStructDir)
                end
                NewName = [CellName, '_', LabelName, ...
                    '_', FileNameSaved, '_Results.mat'];
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