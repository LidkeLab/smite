function [TR, SMD, SMDPreThresh, FileList, TransformList] = batchTrack(obj)
%batchTrack performs batch tracking on multiple sets of SPT data.
% This method is a wrapper around obj.performFullAnalysis() which
% simplifies batch tracking.  That is, this method can be used to simplify
% batch tracking for multiple sets of data.
%
% OUTPUTS:
%   TR: Cell array of Tracking Results structures.
%       (see smi_core.TrackingResults)
%   SMD: Cell array of Single Molecule Data structures.
%        (see smi_core.SingleMoleculeData)
%   SMDPreThresh: Cell array of pre-thresholded 'SMD's.
%   FileList: Cell array of the paths to the files which were used.  The
%             indexing matches that of the cell arrays TR and SMD, e.g.,
%             SMD{n} contains the SMD results from file FileList{n}.
%   TransformList: Cell array of the paths to the transform files used.
%                  The indexing matches that of the cell arrays TR and SMD,
%                  e.g., SMD{n} contains was transformed by
%                  TransformList{n}.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Generate a list of all files in obj.SMF.Data.FileDir which match the 
% pattern given in obj.FilePattern (if needed).
if (obj.FindFiles || isempty(obj.SMF.Data.FileName))
    FileNames = smi_helpers.getFileNames(obj.SMF.Data.FileDir, ...
        obj.FilePattern);
    if isempty(FileNames)
        return
    end
else
    FileNames = obj.SMF.Data.FileName;
end
FileList = fullfile(obj.SMF.Data.FileDir, FileNames);
NFiles = numel(FileNames);
if (obj.Verbose > 1)
    fprintf('smi.SPT.batchTrack(): Found %i files to be tracked.\n', NFiles)
    for ii = 1:NFiles
        fprintf('\t%s\n', FileNames{ii})
    end
end

% Match the files to a channel registration file based on their timestamps.
TransformFiles = ...
    smi_helpers.getFileNames(obj.TransformDir, obj.TransformPattern);
NTransforms = numel(TransformFiles);
if isempty(TransformFiles)
    TransformList = repmat({''}, NFiles, 1);
else
    if (NTransforms > 1)
        TransformTimeStrings = ...
            regexp(TransformFiles, obj.TimeStampRegExp, 'match');
        FileTimeStrings = regexp(FileNames, obj.TimeStampRegExp, 'match');
        PairIndices = smi_helpers.pairTimeStrings(...
            FileTimeStrings, TransformTimeStrings, 'before');
        PairedTransforms = TransformFiles(PairIndices);
        TransformList = fullfile(obj.TransformDir, PairedTransforms);
    else
        % If there's only one transform, we don't need to check timestamps.
        PairedTransforms = repmat(TransformFiles(1), NFiles, 1);
        TransformList = fullfile(obj.TransformDir, PairedTransforms);
    end
    if (obj.Verbose > 1)
        fprintf('smi.SPT.batchTrack(): %i transform(s) matched to data.\n', ...
            NTransforms)
        for ii = 1:NFiles
            fprintf('\t%s <-> %s\n',  FileNames{ii}, PairedTransforms{ii})
        end
    end
end

% Loop through the files and perform tracking. For iterative batch
% tracking, track each file with the current set of SMF parameters before
% making a new estimate from the results.
TR = cell(NFiles, 1);
SMD = cell(NFiles, 1);
SMDPreThresh = cell(NFiles, 1);
IsTestRunInit = obj.IsTestRun;
obj.IsTestRun = (obj.SMF.Tracking.NIterMaxBatch > 1);
IsLastIter = false;
ii = 1;
ParamsHistory = {obj.SMF.Tracking};
while ((ii<=obj.SMF.Tracking.NIterMaxBatch) && ~IsLastIter)
    % Send an update to the command window.
    if (obj.Verbose > 1)
        fprintf(['\tsmi.spt.batchTrack(): ', ...
            'Batch tracking iteration %i...\n'], ii)
    end
    
    % Check parameters for convergence.
    if (ii > 1)
        Converged = obj.checkConvergence(...
            obj.SMF.Tracking, ParamsHistory{ii-1}, ...
            obj.SMF.Tracking.MaxRelativeChange);
        if ((ii==obj.SMF.Tracking.NIterMaxBatch) || Converged)
            % Indicate it's the last iteration and restore the initial
            % setting of obj.IsTestRun (which, if true, will allow for
            % results to be saved on this last iteration).
            obj.IsTestRun = IsTestRunInit;
            IsLastIter = true;
        end
    end
    
    % Track all of the files.
    SMDCat = struct([]);
    for ff = 1:NFiles
        if (obj.Verbose > 0)
            fprintf('smi.SPT.batchTrack(): Tracking file %i of %i...\n', ...
                ff, NFiles)
        end
        obj.SMF.Tracking = ParamsHistory{ii};
        if IsLastIter
            obj.SMF.Tracking.ParamsHistory = ParamsHistory;
        end
        obj.SMF.Data.RegistrationFilePath = TransformList{ff};
        obj.SMF.Data.FileName = FileNames(ff);
        [TR{ff}, SMD{ff}, SMDPreThresh{ff}] = obj.performFullAnalysis();
        SMDCat = smi_core.SingleMoleculeData.catSMD(SMDCat, SMD{ff}, ...
            (obj.Verbose > 1));
    end
    
    % Update the tracking parameters.
    obj.updateTrackingParams(SMDCat)
    ParamsHistory{ii+1, 1} = obj.SMF.Tracking;
    ii = ii + 1;
end
obj.IsTestRun = IsTestRunInit;

% Restore the filenames in SMF (internally, this is overwitten, which
% becomes annoying when testing batch-tracking in the GUI).
obj.SMF.Data.FileName = FileNames;
if (obj.Verbose > 0)
    fprintf('smi.SPT.batchTrack(): Batch-tracking complete.\n')
end


end