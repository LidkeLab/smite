function [SMD, BFStruct] = ...
    driftCorrectBF(SMD, SMF, RefImage, BFStruct, ParamStruct)
%driftCorrectBF performs drift correction from brightfield images.
% This method performs drift correction using brightfield images stored in
% an h5 file.  This is done by estimating shifts between brightfield images
% taken before and after each dataset in SMD, making the drift model a 1D
% polynomial for inter- and intra-dataset drift correction.
%
% INPUTS:
%   SMD: Single Molecule Data structure containing the localizations.
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
%   SMD: Input SMD with applied shifts and updated fields SMD.DriftX and
%        SMD.DriftY.
%   BFStruct: Structure of brightfield images loaded from the data file.

% Created by:
%   David J. Schodt (Lidke Lab 2021)


% Attempt to load the brightfield data (if needed).
if (~exist('BFStruct', 'var') || isempty(BFStruct))
    BFStruct = smi_core.LoadData.readH5File(...
        fullfile(SMF.Data.FileDir, SMF.Data.FileName{1}), 'FocusImages');
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
DefaultParams.ShiftParams = struct();
ParamStruct = smi_helpers.padStruct(ParamStruct, DefaultParams);

% Perform intra-dataset drift correction.
PreSeqImages = zeros([size(RefImage), SMD.NDatasets]);
PostSeqImages = zeros([size(RefImage), SMD.NDatasets]);
for nn = 1:SMD.NDatasets
    PreSeqImages(:, :, nn) = median(BFStruct(nn).Data.PreSeqImages, 3);
    PostSeqImages(:, :, nn) = median(BFStruct(nn).Data.PostSeqImages, 3);
end
SMD = smi_core.DriftCorrection.driftCorrectBFIntra(...
    SMD, PreSeqImages, PostSeqImages, ParamStruct);

% Perform inter-dataset drift correction.
SMD = smi_core.DriftCorrection.driftCorrectBFInter(...
    SMD, RefImage, PreSeqImages, ParamStruct);


end