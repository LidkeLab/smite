function [MLEParams, MLEParamsSE] =  mleOfJumpsBrownian(...
    SquaredDisplacement, FrameLagsAll, LocVarianceSum, NComponents, ...
    FitOptions)
%mleOfJumpsBrownian finds the MLE for Brownian motion jumps.
% This method will estimate the MLE for the likelihood of observed squared
% displacements based on a Brownian motion model with either one diffusing
% population (i.e., one diffusion constant) or two diffusing populations
% (two diffusion constants and the population ratio).
%
% INPUTS:
%   SquaredDisplacement: The squared displacements made by the
%                        trajectory(ies). (NDatax1 numeric array)
%   FrameLagsAll: All of the frame lags associated with the jumps in
%                 'SquaredDisplacement'. (NDatax1 array)
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
%   FitOptions: Fit options sent directly to fminsearch (see doc fminsearch
%               for details)(Default = optimset(@fminsearch)
%               or optimoptions('fmincon') as appropriate)
%
% OUTPUTS:
%   MLEParams: MLE for the desired model parameters..  These will vary
%              based on 'FitMethod'. MLEParams(ii, :) will contain the MLEs
%              corresponding to the data in MSDStruct(ii).
%              (numel(MSDStruct)xNParameters array)
%   MLEParamsSE: Standard errors for the MLEs in 'MLEParams'. These will
%                vary based on 'Method'.
%
% REQUIRES:
%   Optimization Toolbox (for fmincon() when using the N-component models)

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Define default parameters if needed.
if (~exist('NComponents', 'var') || isempty(NComponents))
    NComponents = 2;
end
if (~exist('FitOptions', 'var') || isempty(FitOptions))
    FitOptions = optimoptions('fmincon');
    FitOptions.Display = 'none';
end

% Find the MLE based on the Brownian motion model using a constrained fit.
% NOTE: This model can be found by taking the prob(r|sigma^2=2Dt+loc.error)
%       (which is a product of Gaussians), converting to polar
%       coordinates, and integrating over theta. For multiple frame lags
%       (which we have), we'll also need to integrate over the frame lags
%       times the proportion of each frame lag observed.
% NOTE: The parameter initial guess below was found to be less susceptible
%       to local minima than a "uniform" (each component equal) initial
%       guess (there tends to be a local minimum at the "uniform" initial
%       guess).
CostFunction = @(Params) ...
    -smi_stat.DiffusionEstimator.brownianJumpLikelihood(Params, ...
    SquaredDisplacement, FrameLagsAll, LocVarianceSum);
NFitComponents = 2 * NComponents;
ParamsLowerBound = ...
    [min(SquaredDisplacement./(4*FrameLagsAll))*ones(NComponents, 1);
    zeros(NComponents, 1)];
ParamsUpperBound = ...
    [max(SquaredDisplacement./(4*FrameLagsAll))*ones(NComponents, 1); ...
    ones(NComponents, 1)];
[Prob, Bins] = histcounts(SquaredDisplacement./(4*FrameLagsAll), ...
    NComponents, 'Normalization', 'probability');
DInitGuess = (Bins(1:(end-1))+Bins(2:end)) ./ 2;
AlphaInitGuess = Prob;
ParamsInit = [DInitGuess, AlphaInitGuess].';
Aeq = zeros(NFitComponents);
Aeq(1, (NComponents+1):end) = 1;
beq = zeros(NFitComponents, 1);
beq(1) = 1;
MLEParams = fmincon(CostFunction, ParamsInit, [], [], Aeq, beq, ...
    ParamsLowerBound, ParamsUpperBound, [], FitOptions);

% If requested, estimate the CRLBs of fit parameters.
if (nargout > 1)
    % Compute the CRLB for the variance of the jump distribution.
    r = sqrt(SquaredDisplacement);
    v = 2*FrameLagsAll*MLEParams(1:NComponents).' + mean(LocVarianceSum);
    alpha = MLEParams((NComponents+1):end).';
    Top = alpha .* r .* exp(-0.5*r.^2./v) ...
        .* (0.5*(r.^2)./(v.^3) - (1./v.^2));
    Bottom = sum(alpha .* (r./v) .* exp(-0.5*r.^2./v), 2);
    DTopDV = alpha .* r .* exp(-0.5*r.^2./v) ...
        .* (2./(v.^3) - 2*(r.^2)./(v.^4) + (r.^4)./(4*v.^5));
    DBottomDV = Top;
    DSqDVLogL = sum((Bottom.*DTopDV - Top.*DBottomDV) ./ (Bottom.^2), 1);
    CRLBOfVariance = sqrt(-1 ./ DSqDVLogL);

    % Compute the CRLB for the population ratios.
    DSqDAlphaLogL = sum(-((r./v) .* exp(-0.5*r.^2./v)).^2 ...
        ./ sum(alpha .* (r./v) .* exp(-0.5*r.^2./v), 2).^2, 1);
    CRLBOfAlpha = sqrt(-1 ./ DSqDAlphaLogL);
    MLEParamsSE = [CRLBOfVariance, CRLBOfAlpha].';
end


end