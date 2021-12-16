function [FitParams, FitParamsSE] =  fitCDFOfJumpsBrownian(...
    SortedSquaredDisp, CDFOfJumps, SortedFrameLagsAll, NPoints, ...
    LocVarianceSum, NComponents, Weights, FitMethod, FitOptions)
%fitCDFOfJumpsBrownian fits the CDF of jumps to a Brownian motion model.
% This method will fit the CDF of squared displacements to a Brownian
% motion model with either one diffusing population (i.e., one diffusion
% constant) or two diffusing populations (two diffusion constants and the
% population ratio).
%
% INPUTS:
%   SortedSquaredDisp: The sorted (in ascending order) squared jumps used 
%                      to compute 'CDFOfJumps'. (numeric array)
%   CDFOfJumps: CDF of the jumps (displacements). (NDatax1 array)
%   SortedFrameLagsAll: All of the frame lags associated with the jumps in
%                 'SortedSquaredDisp'. (NDatax1 array)
%   NPoints: The number of data points (or jumps) corresponding to each
%            frame lag in 'FrameLags' (NFrameLagsx1 array)
%   LocVarianceSum: Sum of the localization variances for the two points
%                   used to compute the jumps. This array should be
%                   averaged over x and y.
%                   (NDatax1 numeric array)
%                   NOTE: I don't know which is better: average the
%                         variances, or average the SEs and square them? My
%                         bet is on averaging variances, but I'm not sure.
%                         This can make a big difference in some cases!
%   NComponents: Number of diffusion coefficients to fit.
%                (scalar, integer)(Default = 2)
%   Weights: Weights used for weighted least squares. (NDatax1 array)
%            (Default = ones(NFrameLags, 1) i.e. no weighting)
%   FitMethod: A string specifying the fit method. (Default = 'WeightedLS')
%   FitOptions: Fit options sent directly to fminsearch (see doc fminsearch
%               for details)(Default = optimset(@fminsearch)
%               or optimoptions('fmincon') as appropriate)
%
% OUTPUTS:
%   FitParams: Fit parameters for the MSD fit where
%              MSDFit = FitParams(1) + FitParams(2)*FrameLags
%   FitParamsSE: Standard errors for the MSD fit parameters 'FitParams'.
%
% REQUIRES:
%   Optimization Toolbox (for fmincon() when using the N-component models)

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Define default parameters if needed.
if (~exist('FitMethod', 'var') || isempty(FitMethod))
    FitMethod = 'WeightedLS';
end
NJumps = numel(SortedSquaredDisp);
if (~exist('Weights', 'var') || isempty(Weights))
    Weights = ones(NJumps, 1);
end
if (~exist('NComponents', 'var') || isempty(NComponents))
    NComponents = 2;
end
if (~exist('FitOptions', 'var') || isempty(FitOptions))
    if (NComponents > 1)
        FitOptions = optimoptions('fmincon');
        FitOptions.Display = 'none';
    else
        FitOptions = optimset(@fminsearch);
    end
end

% Fit the CDF of the MSD to the model.
NFitComponents = 2 * NComponents;
switch FitMethod
    case {'LS', 'WeightedLS'}
        % For regular least squares ('LS'), overwrite Weights to just be an
        % array of ones (i.e., no weighting).
        Weights = strcmpi(FitMethod, 'LS')*ones(NJumps, 1) ...
            + strcmpi(FitMethod, 'WeightedLS')*Weights;
        
        % Define lower and upper bounds for the fit parameters (only used
        % for NComponents > 1).
        ParamsLowerBound = ...
            [min(SortedSquaredDisp./(4*SortedFrameLagsAll))*ones(NComponents, 1); ...
            zeros(NComponents, 1)];
        ParamsUpperBound = ...
            [max(SortedSquaredDisp./(4*SortedFrameLagsAll))*ones(NComponents, 1); ...
            ones(NComponents, 1)];
        
        % Define constraints of the form A*x = b (e.g., for now
        % I'm forcing the sum of N population ratios to be == 1).
        Aeq = zeros(NFitComponents);
        Aeq(1, (NComponents+1):end) = 1;
        beq = zeros(NFitComponents, 1);
        beq(1) = 1;
                
        % Fit the CDF of the displacements using least squares. This
        % process is quite slow due to the bootstrap, so I'll only do the
        % bootstrap if the output FitParamsSE was requested.
        % NOTE: This model can be found by taking the
        %       prob(r|sigma^2=2Dt+loc.error) (which is a product of
        %       Gaussians), converting to polar coordinates, integrating
        %       over theta, and then integrating from 0 to r' to get the
        %       CDF. For multiple frame lags (as we have), we'll also need
        %       to integrate over the frame lags times the proportion of
        %       each frame lag observed.
        CostFunction = @(Params, SortedSqJumps, CDFOfJumps) sum(Weights ...
            .* (smi_stat.DiffusionEstimator.brownianJumpCDF(...
            Params, SortedSqJumps, unique(SortedFrameLagsAll), NPoints, ...
            LocVarianceSum) - CDFOfJumps).^2);
        ParamsInit = ...
            [mean(SortedSquaredDisp./(4*SortedFrameLagsAll))*ones(NComponents, 1); ...
            (1/NComponents)*ones(NComponents, 1)];
        if (nargout > 1)
            % For the single component fit, we'll just use fminsearch().
            % For the N-component fit, it sometimes seems important to
            % constrain the fit to get a reasonable result.
            if (NComponents > 1)
                % Perform the constrained fit.
                [FitParams, FitParamsSE] = smi_stat.bootstrapFitCon(...
                    SortedSquaredDisp, CDFOfJumps, ParamsInit, CostFunction, ...
                    [], Aeq, beq, [], [], ...
                    ParamsLowerBound, ParamsUpperBound, FitOptions);
            else
                [FitParams, FitParamsSE] = smi_stat.bootstrapFit(...
                    SortedSquaredDisp, CDFOfJumps, ParamsInit, CostFunction, ...
                    [], FitOptions);
            end
        else
            % For the single component fit, we'll just use fminsearch().
            % For the N-component fit, it sometimes seems important to
            % constrain the fit to get a reasonable result.
            if (NComponents > 1)
                % Perform the constrained fit.
                FitParams = fmincon(@(Params) ...
                    CostFunction(Params, SortedSquaredDisp, CDFOfJumps), ...
                    ParamsInit, [], [], Aeq, beq, ...
                    ParamsLowerBound, ParamsUpperBound);
            else
                FitParams = fminsearch(@(Params) ...
                    CostFunction(Params, SortedSquaredDisp, CDFOfJumps), ...
                    ParamsInit, FitOptions);
            end
        end
    otherwise
        error('Unknown ''FitMethod'' = %s', FitMethod)
end


end