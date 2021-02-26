function [BetaHat, BetaHatSE] = leastSquaresFit(XData, YData, Weights, ...
    CovariateFunctions)
%leastSquaresFit performs a least squares fit on the provided data.
% This function computes a least squares fit for the paired data in XData,
% YData, and Weights.
%
% NOTE: The output 'BetaHatSE' is estimated in the "standard" way presented
%       for linear least squares fitting (at least, it's how I've seen it
%       in all sources I've found). It can be derived by manipulating 
%       var(Beta) = var((X.'*W*X) \ (X.'*W*YData)), estimating var(YData) 
%       by assuming YData(ii) ~ N(Model(ii), Sn(ii)), and then computing 
%       Sn = (1 / (NData-NCovariates)) * (WeightedResiduals).^2.
%
% INPUTS:
%   XData: X data. (NDatax1 array)
%   YData: Y data. (NDatax1 array)
%   Weights: Array of weights used for a weighted least squares fit. These
%            would usually be chosen to be (proportional to) the inverse of
%            the variance of the data in YData.
%            (NDatax1 array)(Default = ones(numel(YData), 1))
%   CovariateFunctions: Cell array of function handles of the "covariates"
%                       (I think covariates might actually mean the results
%                       of applying the function to XData but I'm not
%                       sure). (cell array of function handles)
%                       (Default = {@(X) ones(numel(X), 1); @(X) X},
%                       which is just a line)
%
% OUTPUTS:
%   BetaHat: Estimated fit parameters. (NCovariatesx1 array)
%   BetaHatSE: Estimated standard errors in the fit parameters.
%              (NCovariatesx1 array)

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Set default parameters and reshape arrays as needed.
NData = numel(YData);
if (~exist('Weights', 'var') || isempty(Weights))
    Weights = ones(NData, 1);
end
if (~exist('CovariateFunctions', 'var') || isempty(CovariateFunctions))
    CovariateFunctions = {@(X) ones(numel(X), 1); @(X) X};
end
if isrow(XData)
    XData = XData.';
end
if isrow(YData)
    YData = YData.';
end

% Prepare the X matrix (I'm not sure what it's called, but X is standard
% notation).
NCovariates = numel(CovariateFunctions);
X = zeros(NData, NCovariates);
for ii = 1:NCovariates
    X(:, ii) = CovariateFunctions{ii}(XData);
end

% Perform the least squares fit.
WeightsMatrix = diag(Weights);
XTransposeWXInverse = inv(X.' * WeightsMatrix * X);
BetaHat = XTransposeWXInverse * (X.'*WeightsMatrix*YData);
SumOfWeightedSqResiduals = sum(Weights .* (X*BetaHat-YData).^2);
YDataVarianceEstimate = SumOfWeightedSqResiduals / (NData-NCovariates);
BetaHatSE = sqrt(abs(diag(XTransposeWXInverse * YDataVarianceEstimate)));


end