function [ParamsHat, ParamsHatSE] = bootstrapFit(XData, YData, ...
    ParamsInitialGuess, CostFunction, NBootstrap, FitOptions)
%bootstrapFit minimizes CostFunction and performs a basic boostrap.
% This function takes random samples (with replacement) from XData/YData
% and minimizes the associated cost function to find the parameters from
% that bootstrap. The output  FitParamsSE is the standard deviation of the
% NBootstrap fits. FitParams is the fit obtained for all of the data.
% Note that this approach may not be valid for all cost functions and data.
%
% NOTE: The Statistics and Machine Learning toolbox has several functions
%       for bootstrapping that may be better/more robust than this 
%       function. I've made this function to avoid the need for toolboxes.
%
% INPUTS:
%   XData: X data. (NDatax1 array)
%   YData: Y data. (NDatax1 array)
%   ParamsInitialGuess: Initial guess for the first input Params to the
%                       'CostFunction'. After the first bootstrap sample,
%                       this value will no longer be used (the first result
%                       will take its place).
%   CostFunction: Function handle which takes three inputs Params, XData,
%                 and YData, in that order. Params is an array containing
%                 the parameters that the cost function will be minimized
%                 with respect to. (Function handle)(Default = 
%                 sum((Params(1) + Params(2)*XData - YData).^2))
%   NBootstrap: Number of bootstrap samples to be made from the data.
%               (Default = 100)
%   FitOptions: Fit options sent directly to fminsearch (see doc fminsearch
%               for details)(Default = optimset(@fminsearch))
%
% OUTPUTS:
%   ParamsHat: Estimated fit parameters. (numeric array)
%   ParamsHatSE: Estimated standard errors in the fit parameters.
%                (numeric array)

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Set default parameters and reshape arrays as needed.
if isrow(XData)
    XData = XData.';
end
if isrow(YData)
    YData = YData.';
end
if isrow(ParamsInitialGuess)
    ParamsInitialGuess = ParamsInitialGuess.';
end
if (~exist('CostFunction', 'var') || isempty(CostFunction))
    CostFunction = @(Params, XData, YData) ...
        sum((Params(1) + Params(2)*XData - YData).^2);
end
if (~exist('NBootstrap', 'var') || isempty(NBootstrap))
    NBootstrap = 100;
end
if (~exist('FitOptions', 'var') || isempty(FitOptions))
    FitOptions = optimset(@fminsearch);
end

% Estimate ParamsHat (I'm not sure if this is appropriate. Should I instead
% use ParamsHat = mean(ParamsBootstrap) below? I don't think so, but I
% don't know for sure! -DJS).
CostFunctionSingleInput = @(Params) CostFunction(Params, XData, YData);
ParamsHat = ...
    fminsearch(CostFunctionSingleInput, ParamsInitialGuess, FitOptions);

% Perform the bootstrap.
NData =  numel(YData);
NParams = numel(ParamsInitialGuess);
ParamsBootstrap = zeros(NParams, NBootstrap);
for nn = 1:NBootstrap
    % Select our NData samples (with replacement) from the data and perform
    % the minimization.
    SampleIndices = randi(NData, NData, 1);
    
    % Prepare a cost function that can be passed to fminsearch.
    CostFunctionSingleInput = @(Params) CostFunction(Params, ...
        XData(SampleIndices), YData(SampleIndices));

    % Minimize the cost function.
    ParamsBootstrap(:, nn) = ...
        fminsearch(CostFunctionSingleInput, ParamsHat, FitOptions);
end
ParamsHatSE = std(ParamsBootstrap, [], 2);


end