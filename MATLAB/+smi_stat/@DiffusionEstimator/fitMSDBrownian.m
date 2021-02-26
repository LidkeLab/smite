function [FitParams, FitParamsSE] = ...
    fitMSDBrownian(FrameLags, MSD, NPoints, FitMethod)
%fitMSDBrownian fits an MSD to a Brownian motion model.
% This method will fit a set of mean squared displacement (MSD) to the
% Brownian motion model (i.e., a line).
%
% INPUTS:
%   FrameLags: The frame lags corresponding to the MSD values in 'MSD'.
%              (NFrameLagsx1 array)
%   MSD: The mean squared displacement data corresponding to the frame lags
%        'FrameLags'. (NFrameLagsx1 array)
%   NPoints: The number of displacements used to compute each of MSD.
%            (NFrameLagsx1 array)(only needed if FitMethod = 'WeightedLS')
%            (Default set s.t. weighted least squares weights are all 1)
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
if (~exist('NPoints', 'var') || isempty(NPoints))
    % The weights in the weighted least squares fit will be 
    % NPoints / FrameLags.^2, thus this default makes the weights 1.
    NPoints = FrameLags .^ 2;
end

% Fit the MSD to the model.
switch FitMethod
    case {'LS', 'WeightedLS'}
        % Define the weights to just be 1 (i.e., an unweighted fit) for
        % regular least squares and roughly the inverse of the MSD variance
        % otherwise.
        % NOTE: The weight array is (as I've defined it here)
        %       proportional to the reciprocal of the CRLB of each
        %       point in an MSD plot. The proportionality constant is
        %       dropped because it doesn't affect the weighted fit.
        MSDVariance = NPoints ./ (FrameLags.^2);
        Weights = strcmpi(FitMethod, 'LS')*ones(numel(FrameLags), 1) ...
            + strcmpi(FitMethod, 'WeightedLS')*MSDVariance;
        
        % Perform the least squares fit.
        [BetaHat, BetaHatSE] = ...
            smi_stat.leastSquaresFit(FrameLags, MSD, Weights);
        FitParams = BetaHat.';
        FitParamsSE = BetaHatSE.';
    otherwise
        error('Unknown ''FitMethod'' = %s', FitMethod)
end


end