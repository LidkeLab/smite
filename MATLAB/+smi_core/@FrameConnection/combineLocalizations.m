function [SMDCombined] = combineLocalizations(SMD, SMF)
%combineLocalizations combines localizations in SMD with the same ConnectID
% This method combines localizations in SMD which share the same value of
% SMD.ConnectID.  That is, this method combines the frame-connected
% localizations into a single localization with higher precision.
%
% INPUTS:
%   SMD: Single Molecule Data structure with a populated SMD.ConnectID.
%   SMF: Single Molecule Fitting structure (needed for the field
%        SMF.Fitting.FitType).  If not provided, a best guess for the fit
%        type will be inferred from the data.
%
% OUTPUTS:
%   SMDCombined: Single Molecule Data structure with combined
%                localizations, where SMD.ConnectID is now unique for every
%                localization.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Isolate/organize some SMD arrays.
ConnectID = SMD.ConnectID;
[ConnectID, SortIndices] = sort(ConnectID);
X = SMD.X(SortIndices);
Y = SMD.Y(SortIndices);
X_SE = SMD.X_SE(SortIndices);
Y_SE = SMD.Y_SE(SortIndices);
FrameNum = SMD.FrameNum(SortIndices);
DatasetNum = SMD.DatasetNum(SortIndices);
Photons = SMD.Photons(SortIndices);
Bg = SMD.Bg(SortIndices);
LogLikelihood = SMD.LogLikelihood(SortIndices);
NLocalizations = numel(FrameNum);
if (~exist('SMF', 'var') || isempty(SMF))
    % Try to determine the fit type based on the presence and size of some
    % array fields.
    if (isfield(SMD, 'PSFSigma_SE') ...
            && (numel(SMD.PSFSigma_SE)==NLocalizations))
        FitType = 'XYNBS';
        PSFSigma = SMD.PSFSigma(SortIndices);
        PSFSigma_SE = SMD.PSFSigma_SE(SortIndices);
    elseif (isfield(SMD, 'PSFSigmaX_SE') ...
            && (numel(SMD.PSFSigmaX_SE)==NLocalizations))
        FitType = 'XYNBSXSY';
        PSFSigmaX = SMD.PSFSigmaX(SortIndices);
        PSFSigmaY = SMD.PSFSigmaY(SortIndices);
        PSFSigmaX_SE = SMD.PSFSigmaX_SE(SortIndices);
        PSFSigmaY_SE = SMD.PSFSigmaY_SE(SortIndices);
    elseif (isfield(SMD, 'Z_SE') && (numel(SMD.Z_SE)==NLocalizations))
        FitType = 'XYZNB';
        Z = SMD.Z(SortIndices);
        Z_SE = SMD.Z_SE(SortIndices);
    else
        FitType = 'XYNB';
    end
else
    FitType = SMF.Fitting.FitType;
end

% Loop over the unique connect IDs and combine the associated
% localizations.
NLocPerID = groupcounts(ConnectID);
NLocCumulative = [0; cumsum(NLocPerID)];
UniqueIDs = unique(ConnectID);
NUnique = numel(UniqueIDs);
SMDCombined = smi_core.SingleMoleculeData.isolateSubSMD(SMD, 1);
for ii = 1:NUnique
    IndexArray = (1:NLocPerID(ii)).' + NLocCumulative(ii);
    SMDCombined.X(ii, 1) = sum(X(IndexArray)./X_SE(IndexArray).^2) ...
        / sum(1./X_SE(IndexArray).^2);
    SMDCombined.Y(ii, 1) = sum(Y(IndexArray)./Y_SE(IndexArray).^2) ...
        / sum(1./Y_SE(IndexArray).^2);
    SMDCombined.X_SE(ii, 1) = sqrt(1 ./ sum(1./X_SE(IndexArray).^2));
    SMDCombined.Y_SE(ii, 1) = sqrt(1 ./ sum(1./Y_SE(IndexArray).^2));
    SMDCombined.ConnectID(ii, 1) = ConnectID(IndexArray(1));
    SMDCombined.FrameNum(ii, 1) = FrameNum(IndexArray(NLocPerID(ii)));
    SMDCombined.DatasetNum(ii, 1) = DatasetNum(IndexArray(1));
    SMDCombined.Photons(ii, 1) = sum(Photons(IndexArray));
    SMDCombined.Bg(ii, 1) = sum(Bg(IndexArray));
    SMDCombined.LogLikelihood(ii, 1) = sum(LogLikelihood(IndexArray));
end
SMDCombined.NCombined(:, 1) = NLocPerID;

% If the fit type wasn't 'XYNB', we still need to combine some other fields
% (I'm doing this in a separate loop down here for speed purposes, since we
% usually won't need to do this).
switch FitType
    case 'XYNB'
        return
    case 'XYNBS'
        for ii = 1:NUnique
            IndexArray = (1:NLocPerID(ii)).' + NLocCumulative(ii);
            SMDCombined.PSFSigma(ii, 1) = ...
                sum(PSFSigma(IndexArray)./PSFSigma_SE(IndexArray).^2) ...
                / sum(1./PSFSigma_SE(IndexArray).^2);
            SMDCombined.PSFSigma_SE(ii, 1) = ...
                sqrt(1 / sum(1./PSFSigma_SE(IndexArray).^2));
        end
    case 'XYNBSXSY'
        for ii = 1:NUnique
            IndexArray = (1:NLocPerID(ii)).' + NLocCumulative(ii);
            SMDCombined.PSFSigmaX(ii, 1) = ...
                sum(PSFSigmaX(IndexArray)./PSFSigmaX_SE(IndexArray).^2) ...
                / sum(1./PSFSigmaX_SE(IndexArray).^2);
            SMDCombined.PSFSigmaY(ii, 1) = ...
                sum(PSFSigmaY(IndexArray)./PSFSigmaY_SE(IndexArray).^2) ...
                / sum(1./PSFSigmaY_SE(IndexArray).^2);
            SMDCombined.PSFSigmaX_SE(ii, 1) = ...
                sqrt(1 / sum(1./PSFSigmaX_SE(IndexArray).^2));
            SMDCombined.PSFSigmaY_SE(ii, 1) = ...
                sqrt(1 / sum(1./PSFSigmaY_SE(IndexArray).^2));
        end
    case 'XYZNB'
        for ii = 1:NUnique
            IndexArray = (1:NLocPerID(ii)).' + NLocCumulative(ii);
            SMDCombined.Z(ii, 1) = ...
                sum(Z(IndexArray)./Z_SE(IndexArray).^2) ...
                / sum(1./Z_SE(IndexArray).^2);
            SMDCombined.Z_SE(ii, 1) = ...
                sqrt(1 / sum(1./Z_SE(IndexArray).^2));
        end
end


end