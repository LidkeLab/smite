function [FitParams, FitParamsSE] =  fitCDFOfJumpsBrownian(...
    SortedJumps, CDFOfJumps, FrameLags, NPoints, ...
    Weights, FitMethod)
%fitCDFOfJumpsBrownian fits the CDF of jumps to a Brownian motion model.
% This method will fit the CDF of an MSD to the Brownian motion model 
% (i.e., a line).
%
% INPUTS:
%   SortedJumps: The sorted jumps values used to compute 'CDFOfJumps' in 
%                ascending order. (NDatax1 numeric array)
%   CDFOfJumps: CDF of the jumps (displacements). (NDatax1 array)
%   FrameLags: All of the unique frame lags associated with the jumps in
%              'SortedJumps'. Note that this isn't necessarily the same
%              size as SortedJumps, since SortedJumps can contain multiple
%              jumps for each frame lag.
%              (NFrameLagsx1 array)
%   NPoints: The number of data points (or jumps) corresponding to each 
%            frame lag in 'FrameLags' (NFrameLagsx1 array)
%   Weights: Weights used for weighted least squares. (NDatax1 array)
%            (Default = ones(NFrameLags, 1) i.e. no weighting)
%   FitMethod: A string specifying the fit method. (Default = 'LS')
%
% OUTPUTS:
%   FitParams: Fit parameters for the MSD fit where
%              MSDFit = FitParams(1) + FitParams(2)*FrameLags
%   FitParamsSE: Standard errors for the MSD fit parameters 'FitParams'.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Define default parameters if needed.
if (~exist('FitMethod', 'var') || isempty(FitMethod))
    FitMethod = 'LS';
end
NJumps = numel(SortedJumps);
if (~exist('Weights', 'var') || isempty(Weights))
    Weights = ones(NJumps, 1);
end

% Fit the CDF of the MSD to the model.
switch FitMethod
    case {'LS', 'WeightedLS'}
        % For regular least squares ('LS'), overwrite Weights to just be an
        % array of ones (i.e., no weighting).
        Weights = strcmpi(FitMethod, 'LS')*ones(NJumps, 1) ...
            + strcmpi(FitMethod, 'WeightedLS')*Weights;
        
        % Fit the CDF of the displacements using least squares. This
        % process is quite slow due to the bootstrap, so I'll only do the
        % bootstrap if the output FitParamsSE was requested.
        % NOTE: This model can be found by taking the prob(r|sigma^2=2Dt)
        %       (which is a product of Gaussians), converting to polar
        %       coordinates, integrating over theta, and then integrating
        %       from 0 to r' to get the CDF. For multiple frame lags (as we
        %       have), we'll also need to integrate over the frame lags
        %       times the proportion of each frame lag observed.
        CostFunction = @(Params, SortedJumps, CDFOfJumps) sum(Weights ...
            .* (smi_stat.DiffusionEstimator.brownianJumpCDF(...
            Params, SortedJumps, FrameLags, NPoints) ...
            - CDFOfJumps).^2);
        if (nargout > 1)
            [FitParams, FitParamsSE] = smi_stat.bootstrapFit(...
                SortedJumps, CDFOfJumps, 0.1, CostFunction);
        else
            FitParams = fminsearch(@(Params) ...
                CostFunction(Params, SortedJumps, CDFOfJumps), 0.1);
        end
    otherwise
        error('Unknown ''FitMethod'' = %s', FitMethod)
end


end