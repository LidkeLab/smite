function [FitParams, FitParamsSE] = ...
    fitCDFOfJumpsBrownian(SortedJumps, FrameLags, CDFOfJumps, FitMethod)
%fitCDFOfJumpsBrownian fits the CDF of jumps to a Brownian motion model.
% This method will fit the CDF of an MSD to the Brownian motion model 
% (i.e., a line).
%
% INPUTS:
%   SortedJumps: The sorted jumps values used to compute 'CDFOfJumps' in 
%                ascending order. (NFrameLagsx1 array)
%   FrameLags: Frame lags corresponding to the jumps in 'SortedJumps'.
%              (NFrameLagsx1 array)
%   CDFOfJumps: CDF of the jumps (displacements). (NFrameLagsx1 array)
%   FitMethod: A string specifying the fit method. (Default = 'WeightedLS')
%
% OUTPUTS:
%   FitParams: Fit parameters for the MSD fit where
%              MSDFit = FitParams(1) + FitParams(2)*FrameLags
%   FitParamsSE: Standard errors for the MSD fit parameters 'FitParams'.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Define default parameters if needed.
if (~exist('FitMethod', 'var') || isempty(FitMethod))
    FitMethod = 'WeightedLS';
end

% Fit the CDF of the MSD to the model.
switch FitMethod
    case 'LS'
        % Fit the CDF of the displacements using least squares.
        % NOTE: In this model, X is 2*NDimensions*DiffusionConstant
        ModelFunction = @(X) ...
            1 - exp(-SortedJumps ./ (X*FrameLags));
        CostFunction = @(X) ...
            sum((CDFOfJumps - ModelFunction(X)).^2);
        FitParams = fminsearch(CostFunction, 0.1);
        FitParamsSE = NaN;
%         [BetaHat, BetaHatSE] = ...
%             smi_stat.leastSquaresFit(FrameLags, MSD);
%         FitParams = BetaHat.';
%         FitParamsSE = BetaHatSE.';
    case 'WeightedLS'
%         % Fit the MSD using weighted linear least squares.
%         % NOTE: The weight array is (as I've defined it here)
%         %       proportional to the reciprocal of the CRLB of each
%         %       point in an MSD plot. The proportionality constant is
%         %       dropped because it doesn't affect the weighted fit.
%         Weights = NPoints ./ (FrameLags.^2);
%         [BetaHat, BetaHatSE] = ...
%             smi_stat.leastSquaresFit(FrameLags, MSD, Weights);
%         FitParams = BetaHat.';
%         FitParamsSE = BetaHatSE.';
    otherwise
        error('Unknown ''FitMethod'' = %s', FitMethod)
end


end