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


% Generate a list of all files in obj.SMF.FileDir which match the pattern
% given in obj.FilePattern.
FileNames = smi_helpers.getFileNames(obj.SMF.Data.FileDir, obj.FilePattern);
NFiles = numel(FileNames);
TR = cell(NFiles, 1);
SMD = cell(NFiles, 1);
SMDPreThresh = cell(NFiles, 1);
FileList = cell(NFiles, 1);
TransformList = cell(NFiles, 1);
if isempty(FileNames)
    return
end
FileList = fullfile(obj.SMF.Data.FileDir, FileNames);
if (obj.Verbose > 1)
    fprintf('smi.SPT.batchTrack(): Found %i files to be tracked.\n', NFiles)
    for ii = 1:NFiles
        fprintf('\t%c\n', FileNames{ii})
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
        PairedTransforms = repmat(TransformFiles{1}, NFiles, 1);
        TransformList = fullfile(obj.TransformDir, PairedTransforms);
    end
    if (obj.Verbose > 1)
        fprintf('smi.SPT.batchTrack(): %i transforms matched to files.\n', ...
            NTransforms)
        for ii = 1:NFiles
            fprintf('\t%c <-> %c\n',  FileNames{ii}, PairedTransforms{ii})
        end
    end
end

% Loop through the files and perform tracking.
for ff = 1:NFiles
    if (obj.Verbose > 0)
        fprintf('smi.SPT.batchTrack(): tracking file %i of %i.\n', ...
            ff, NFiles)
    end
    obj.SMF.Data.RegistrationFilePath = TransformList{ff};
    obj.SMF.Data.FileName = FileNames(ff);
    [TR{ff}, SMD{ff}, SMDPreThresh{ff}] = obj.performFullAnalysis();
end
if (obj.Verbose > 0)
    fprintf('smi.SPT.batchTrack(): batch-tracking complete.\n')
end


end