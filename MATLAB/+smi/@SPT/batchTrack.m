function [TR, SMD, FileList, TransformList] = batchTrack(obj)
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
FileList = cell(NFiles, 1);
TransformList = cell(NFiles, 1);
if isempty(FileNames)
    return
end
FileList = fullfile(obj.SMF.Data.FileDir, FileNames);

% Generate a list of all channel registration files in obj.TransformDir
% which match the pattern given in obj.TransformPattern.  Also, place their
% timestamps in an array for later use.
TransformFiles = ...
    smi_helpers.getFileNames(obj.TransformDir, obj.TransformPattern);
NTransforms = numel(TransformFiles);
TransformTimes = zeros(NTransforms, 1);
for ff = 1:NTransforms
    TimeStamp = regexp(TransformFiles{ff}, obj.TimeStampRegExp, 'match');
    TransformTimes(ff) = smi_helpers.convertTimeStringToNum(...
        TimeStamp{1}, obj.TimeStampDelimiter);
end

% Create the TransformList by matching timestamps in FileNames to the
% nearest timestamps in TransformFiles.
TransformList = cell(NFiles, 1);
if isempty(TransformFiles)
    TransformList = repmat({''}, NFiles, 1);
else
    for ff = 1:NFiles
        % Find the timestamp in the filename.
        TimeStamp = regexp(FileNames{ff}, obj.TimeStampRegExp, 'match');
        TimeNum = smi_helpers.convertTimeStringToNum(...
            TimeStamp{1}, obj.TimeStampDelimiter);
        
        % Compare the timestamp of the file to those in the transform files
        % and select the closest match.
        % NOTE: As written, this will find the transform timestamped
        %       closest in time BEFORE the data was taken.
        TimeDiff = TimeNum - TransformTimes;
        ValidInd = find(TimeDiff >= 0);
        [~, MinDiffInd] = min(TimeDiff(ValidInd));
        TransformList{ff} = fullfile(obj.TransformDir, ...
            TransformFiles{ValidInd(MinDiffInd)});
    end
end

% Loop through the files and perform tracking.
ResultsDir = obj.SMF.Data.ResultsDir;
for ff = 1:NFiles
    obj.SMF.Data.RegistrationFilePath = TransformList{ff};
    obj.SMF.Data.FileName = FileNames(ff);
    obj.SMF.Data.ResultsDir = ResultsDir;
    [TR{ff}, SMD{ff}] = obj.performFullAnalysis();
end


end