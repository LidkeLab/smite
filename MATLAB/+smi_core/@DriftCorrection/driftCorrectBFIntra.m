function [SMD] = ...
    driftCorrectBFIntra(SMD, PreSeqImages, PostSeqImages, ParamStruct)
%driftCorrectBFIntra computes intra-DS drift correction from brightfield.
% This method computes intra-dataset drift correction using brightfield 
% images stored in an h5 file.
%
% INPUTS:
%   SMD: Single Molecule Data structure containing the localizations.
%   PreSeqImages: Stack of brightfield images collected before each
%                 dataset.
%   PostSeqImages: Stack of brightfield images collected after each
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

% Perform intra-dataset drift correction.
if (~isfield(SMD, 'DriftX') || ~isfield(SMD, 'DriftY') ...
        || isempty(SMD.DriftX) || isempty(SMD.DriftY))
    SMD.DriftY = NaN(SMD.NFrames, SMD.NDatasets);
    SMD.DriftX = NaN(SMD.NFrames, SMD.NDatasets);
end
for nn = 1:size(PreSeqImages, 3)
    % Compute the shift between the pre- and post-sequence images.
    Shift = smi_stat.findOffsetIter(...
        PreSeqImages(:, :, nn), PostSeqImages(:, :, nn), ...
        ParamStruct.NIterMax, ParamStruct.Tolerance, ...
        ParamStruct.CorrParams, ParamStruct.ShiftParams);

    % Apply and store the shift in SMD.
    CurrentDS = (SMD.DatasetNum == nn);
    FrameNum = SMD.FrameNum(CurrentDS);
    SMD.Y(CurrentDS) = SMD.Y(CurrentDS) + (Shift(1)/SMD.NFrames)*(FrameNum-1);
    SMD.X(CurrentDS) = SMD.X(CurrentDS) + (Shift(2)/SMD.NFrames)*(FrameNum-1);
    SMD.DriftY(:, nn) = (Shift(1)/SMD.NFrames)*(0:(SMD.NFrames-1)).';
    SMD.DriftX(:, nn) = (Shift(2)/SMD.NFrames)*(0:(SMD.NFrames-1)).';
end


end