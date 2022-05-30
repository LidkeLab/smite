function [SMDCombined] = combineLocalizations(SMD, SMF)
%combineLocalizations combines localizations in SMD with the same ConnectID
% This method combines localizations in SMD which share the same value of
% SMD.ConnectID.  That is, this method combines the frame-connected
% localizations into a single localization with higher precision.
%
% INPUTS:
%   SMD: Single Molecule Data structure with a populated SMD.ConnectID.
%   SMF: Single Molecule Fitting structure (needed for the field
%        SMF.Fitting.FitType).
%
% OUTPUTS:
%   SMDCombined: Single Molecule Data structure with combined
%                localizations, where SMD.ConnectID is now unique for every
%                localization.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Define a default for 'FitType' if needed.
if (~exist('SMF', 'var') || isempty(SMF))
    SMF = smi_core.SingleMoleculeFitting;
end

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
switch SMF.Fitting.FitType
    case 'XYNBS'
        PSFSigma = SMD.PSFSigma(SortIndices);
        PSFSigma_SE = SMD.PSFSigma_SE(SortIndices);
    case 'XYNBSXSY'
        PSFSigmaX = SMD.PSFSigmaX(SortIndices);
        PSFSigmaY = SMD.PSFSigmaY(SortIndices);
        PSFSigmaX_SE = SMD.PSFSigmaX_SE(SortIndices);
        PSFSigmaY_SE = SMD.PSFSigmaY_SE(SortIndices);
    case 'XYZNB'
        Z = SMD.Z(SortIndices);
        Z_SE = SMD.Z_SE(SortIndices(1:numel(SMD.Z_SE)));
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
end
SMDCombined.NCombined = NLocPerID;

% If the fit type wasn't 'XYNB', we still need to combine some other fields
% (I'm doing this in a separate loop down here for speed purposes, since we
% usually won't need to do this).
switch SMF.Fitting.FitType
    case 'XYNB'
        % For fit type 'XYNB', every value of PSFSigma should be the same.
        SMDCombined.PSFSigma = ...
            SMD.PSFSigma(1) * ones(size(SMDCombined.FrameNum));
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

% Combine some fields that may or may not be in the SMD (e.g., a simulation
% might not produce a value for log-likelihood, which can become a burden
% to deal with).
% NOTE: For Photons_SE and Bg_SE, I'm not sure if this is the correct
%       result.  However, it's probably decent: these are the correct
%       results assuming only Poisson (shot) noise.
if ~isempty(SMD.Photons_SE)
    Photons_SE = SMD.Photons_SE(SortIndices);
    for ii = 1:NUnique
        IndexArray = (1:NLocPerID(ii)).' + NLocCumulative(ii);
        SMDCombined.Photons_SE(ii, 1) = sqrt(sum(Photons_SE(IndexArray).^2));
    end
end
if ~isempty(SMD.Bg_SE)
    Bg_SE = SMD.Bg_SE(SortIndices);
    for ii = 1:NUnique
        IndexArray = (1:NLocPerID(ii)).' + NLocCumulative(ii);
        SMDCombined.Bg_SE(ii, 1) = sqrt(sum(Bg_SE(IndexArray).^2));
    end
end
if ~isempty(SMD.LogLikelihood)
    LogLikelihood = SMD.LogLikelihood(SortIndices);
    for ii = 1:NUnique
        IndexArray = (1:NLocPerID(ii)).' + NLocCumulative(ii);
        SMDCombined.LogLikelihood(ii, 1) = sum(LogLikelihood(IndexArray));
    end
end
if ~isempty(SMD.PValue)
    SMDCombined.PValue = smi_core.GaussMLE.pValue(SMF.Fitting.NParams, ...
        SMF.BoxFinding.BoxSize, ...
        SMDCombined.LogLikelihood);
end
if ~isempty(SMD.ThreshFlag)
    ThreshFlag = SMD.ThreshFlag(SortIndices);
    for ii = 1:NUnique
        IndexArray = (1:NLocPerID(ii)).' + NLocCumulative(ii);
        SMDCombined.ThreshFlag(ii, 1) = sum(ThreshFlag(IndexArray));
    end
end
if ~(isempty(SMD.XBoxCorner) || isempty(SMD.YBoxCorner))
    XBoxCorner = SMD.XBoxCorner(SortIndices);
    YBoxCorner = SMD.YBoxCorner(SortIndices);
    for ii = 1:NUnique
        SMDCombined.XBoxCorner(ii, 1) = ...
            XBoxCorner(NLocPerID(ii) + NLocCumulative(ii), 1);
        SMDCombined.YBoxCorner(ii, 1) = ...
            YBoxCorner(NLocPerID(ii) + NLocCumulative(ii), 1);
    end
end


end