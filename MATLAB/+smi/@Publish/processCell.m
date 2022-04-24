function processCell(obj, CellName)
%processCell will process data corresponding to CellName.
% This method will find the sub-directories for the cell CellName, which
% themselves contain the data for each label of the acquisition, and
% analyze that data.
%
% INPUTS:
%   CellName: Char. array/string of the Cell* directory (i.e., the distinct
%             Cell* names in the directory structure
%             *\Cell*\Label*\Data*.h5)

% Created by:
%   David J. Schodt (Lidke Lab 2021)


% Determine the names of the sub-directories of interest within
% CellName.  These correspond to single labels imaged during the
% super-resolution experiment.
LabelNames = smi_helpers.getDirectoryNames(...
    fullfile(obj.CoverslipDir, CellName), 'Label*');
NLabels = numel(LabelNames);
if (obj.Verbose > 1)
    fprintf('\tPublish.processCell(): %i label directories found:\n', ...
        NLabels)
    for ii = 1:NLabels
        fprintf('\t\t%s\n', LabelNames{ii})
    end
end

% Loop through each of the label directories and process the data.  If the
% processing fails on a given label ii, proceed with the next label anyways
% (these results might still be useful).
obj.SMLM = smi.SMLM(copy(obj.SMF));
obj.SMLM.Verbose = obj.Verbose;
obj.SMLM.SRImageZoom = obj.SRImageZoom;
obj.SMLM.SRCircImZoom = obj.SRCircleImageZoom;
for ii = 1:NLabels
    % If LabelID was specified, skip all labels except those which exist in
    % LabelID.  However, if obj.LabelID is empty, then we wish to analyze
    % all LabelID's available.
    if ~(ismember(ii, obj.LabelID) || isempty(obj.LabelID))
        continue
    end

    % Attempt to process the data for label ii.
    try
        obj.processLabel(CellName, LabelNames{ii});
    catch MException
        if obj.Verbose
            warning(['Publish.processCell(): ', ...
                'Processing failed for %s\n%s, %s'], ...
                fullfile(CellName, LabelNames{ii}), ...
                MException.identifier, MException.message)
        end

        % Store the error information in the log file.
        obj.ErrorLog = [obj.ErrorLog; ...
            {CellName, LabelNames{ii}, MException}];
        ErrorLog = obj.ErrorLog;
        save(obj.LogFilePath, 'ErrorLog', '-append')
    end
end

% If all labels for this cell were processed successfully, create an
% overlay image of the multiple labels, storing the overlay in the top
% level directory for easy access.
if (obj.GenerateSR && (NLabels>1))
    % Prepare overlay masks.
    CellNumber = regexp(CellName, '\d*', 'match');
    CellNumber = str2double(CellNumber{1});
    Mask = obj.genBFMask(obj.FocusImageStructs(CellNumber, :), ...
        obj.MaxBrightfieldShift);
    MaskName = sprintf('%inm', ...
        round(obj.MaxBrightfieldShift * obj.SMF.Data.PixelSize * 1e3));
    if (obj.MaxBrightfieldShift > 0)
        save(fullfile(obj.SaveBaseDir, ...
            sprintf('%s_%s_Mask.mat', CellName, MaskName)), 'Mask')
    end

    % Prepare the overlays.
    try
        obj.genSROverlays(...
            fullfile(obj.SaveBaseDir, CellName), ...
            obj.SaveBaseDir, obj.SMF.Data.AnalysisID, ...
            Mask, MaskName);
    catch MException
        if obj.Verbose
            warning(['Publish.processCell(): Overlay image ', ...
                'generation failed for %s\n%s, %s'], ...
                CellName, MException.identifier, MException.message)
        end

        % Store the error information in the log file.
        obj.ErrorLog = [obj.ErrorLog; ...
            {CellName, LabelNames{1}, MException}];
        ErrorLog = obj.ErrorLog;
        save(obj.LogFilePath, 'ErrorLog', '-append')
    end
end


end