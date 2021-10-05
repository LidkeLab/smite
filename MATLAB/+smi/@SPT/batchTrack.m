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
        TransformList = ...
            fullfile(obj.TransformDir, TransformFiles(PairIndices));
    else
        % If there's only one transform, we don't need to check timestamps.
        TransformList = repmat(...
            {fullfile(obj.TransformDir, TransformFiles{1})}, ...
            NFiles, 1);
    end    
end

% Loop through the files and perform tracking.
for ff = 1:NFiles
    obj.SMF.Data.RegistrationFilePath = TransformList{ff};
    obj.SMF.Data.FileName = FileNames(ff);
    [TR{ff}, SMD{ff}, SMDPreThresh{ff}] = obj.performFullAnalysis();
end


end