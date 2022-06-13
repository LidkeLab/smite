function [RefImage, PreSeqImages, PostSeqImages, ParamStruct] = ...
    driftCorrectBFInit(NDatasets, SMF, RefImage, BFStruct, ParamStruct)
%driftCorrectBF initializes drift correction from brightfield images
% This method performs drift correction using brightfield images stored in
% an h5 file.  This is done by estimating shifts between brightfield images
% taken before and after each dataset in SMD, making the drift model a 1D
% polynomial for inter- and intra-dataset drift correction.
%
% INPUTS:
%   SMD: Single Molecule Data structure containing the localizations.
%   NDatasets: Total number of datasets to be processed.
%   SMF: Single Molecule Fitting structure defining the path to the raw
%        data file.
%   RefImage: Reference image to which all datasets will be corrected to.
%             (Default = pre dataset 1 focus image).
%   BFStruct: Structure of brightfield images (see default setting below).
%   ParamStruct: Structure of parameters sent to smi_stat.findOffsetIter().
%                (fields NIterMax, Tolerance, CorrParams, ShiftParams are
%                passed as inputs to smi_stat.findOffsetIter())
%
% OUTPUTS:
%   RefImage: Reference image to which all datasets will be corrected to.
%   PreSeqImages: 
%   PostSeqImages: 
%   ParamStruct: Structure of parameters sent to smi_stat.findOffsetIter().

% Created by:
%   David J. Schodt (Lidke Lab 2021) and Michael J. Wester (2022)

% Attempt to load the brightfield data (if needed).
if (~exist('BFStruct', 'var') || isempty(BFStruct))
    try
        FilePath = fullfile(SMF.Data.FileDir, SMF.Data.FileName{1});
        H5FileStruct = h5info(FilePath);
        FileGroupList = {H5FileStruct.Groups.Groups.Groups(1).Groups.Name};
        FocusImagesPresent = any(contains(FileGroupList, 'FocusImages'));
        if FocusImagesPresent
           BFStruct = smi_core.LoadData.readH5File(FilePath, 'FocusImages');
        else
           error(sprintf('Cannot extract group ''FocusImages'' from %s', ...
                 FilePath));
        end
    catch ME
        error(sprintf('Cannot extract group ''FocusImages'' from %s', ...
              FilePath));
    end
end

% Set a defaults if needed.
if (~exist('RefImage', 'var') || isempty(RefImage))
    RefImage = median(BFStruct(1).Data.PreSeqImages, 3);
end
if (~exist('ParamStruct', 'var') || isempty(ParamStruct))
    ParamStruct = struct();
end
DefaultParams.NIterMax = 3;
DefaultParams.Tolerance = [];
DefaultParams.CorrParams = struct();
DefaultParams.CorrParams.SuppressWarnings = true;
DefaultParams.ShiftParams = struct();
ParamStruct = smi_helpers.padStruct(ParamStruct, DefaultParams);

% Perform intra-dataset drift correction.
PreSeqImages = zeros([size(RefImage), NDatasets]);
PostSeqImages = zeros([size(RefImage), NDatasets]);
for nn = 1:NDatasets
    PreSeqImages(:, :, nn) = median(BFStruct(nn).Data.PreSeqImages, 3);
    PostSeqImages(:, :, nn) = median(BFStruct(nn).Data.PostSeqImages, 3);
end

% Moved to SMLM.
% Perform intra-dataset drift correction.
%SMD = smi_core.DriftCorrection.driftCorrectBFIntra(...
%    SMD, PreSeqImages, PostSeqImages, ParamStruct);
%
% Perform inter-dataset drift correction.
%SMD = smi_core.DriftCorrection.driftCorrectBFInter(...
%    SMD, RefImage, PreSeqImages, ParamStruct);


end
