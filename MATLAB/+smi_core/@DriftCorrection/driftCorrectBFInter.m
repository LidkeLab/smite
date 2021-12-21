function [SMD] = ...
    driftCorrectBFInter(SMD, RefImage, PreSeqImages, ParamStruct)
%driftCorrectBFInter computes inter-DS drift correction from brightfield.
% This method computes inter-dataset drift correction using brightfield 
% images stored in an h5 file.
%
% INPUTS:
%   SMD: Single Molecule Data structure containing the localizations.
%   RefImage: Reference image to which all datasets will be corrected to.
%             (Default = pre dataset 1 focus image).
%   PreSeqImages: Stack of brightfield images collected before each
%                 dataset.
%   ParamStruct: Structure of parameters sent to smi_stat.findOffsetIter().
%                (fields NIterMax, Tolerance, CorrParams, ShiftParams are
%                passed as inputs to smi_stat.findOffsetIter())
%
% OUTPUTS:
%   SMD: Input SMD with applied shifts and updated fields SMD.DriftX and
%        SMD.DriftY.

% Created by:
%   David J. Schodt (Lidke Lab 2021)


% Set a defaults if needed.
if (~exist('ParamStruct', 'var') || isempty(ParamStruct))
    ParamStruct = struct();
end
DefaultParams.NIterMax = 3;
DefaultParams.Tolerance = [];
DefaultParams.CorrParams = struct();
DefaultParams.ShiftParams = struct();
ParamStruct = smi_helpers.padStruct(ParamStruct, DefaultParams);

% Perform inter-dataset drift correction.
if (~isfield(SMD, 'DriftX') || ~isfield(SMD, 'DriftY') ...
        || isempty(SMD.DriftX) || isempty(SMD.DriftY))
    SMD.DriftY = zeros(SMD.NFrames, SMD.NDatasets);
    SMD.DriftX = zeros(SMD.NFrames, SMD.NDatasets);
end
for nn = 1:size(PreSeqImages, 3)
    % Compute the shift between the pre-sequence images and the reference.
    Shift = smi_stat.findOffsetIter(...
        RefImage, PreSeqImages(:, :, nn), ...
        ParamStruct.NIterMax, ParamStruct.Tolerance, ...
        ParamStruct.CorrParams, ParamStruct.ShiftParams);

    % Apply and store the shift in SMD.
    CurrentDS = (SMD.DatasetNum == nn);
    SMD.Y(CurrentDS) = SMD.Y(CurrentDS) + Shift(1);
    SMD.X(CurrentDS) = SMD.X(CurrentDS) + Shift(2);
    SMD.DriftY(:, nn) = SMD.DriftY(:, nn) + Shift(1);
    SMD.DriftX(:, nn) = SMD.DriftX(:, nn) + Shift(2);
end


end