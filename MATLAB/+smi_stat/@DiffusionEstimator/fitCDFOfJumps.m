function [FitParams, FitParamsSE] = fitCDFOfJumps(MSDStruct, FitMethod, ...
    NComponents, DiffusionModel, Verbose)
%fitCDFOfJumps fits the CDF of displacement data (jumps).
% This method will fit a model to the CDF (CPD) of trajectory
% displacements.
%
% INPUTS:
%   MSDStruct: A structure array of MSD data as output from computeMSD()
%              (see computeMSD() for more details). The MSDStruct must
%              already be populated with SortedJumps and CDFOfJumps (see
%              computeCDFOfMSD()).
%   FitMethod: A string specifying the fit method. (Default = 'WeightedLS')
%   NComponents: Number of diffusive components used in fit. (Default = 1)
%   DiffusionModel: A string specifying the diffusion model to fit to the
%                   MSD. See options in DiffusionEstimator class property
%                   'DiffusionModel'. (Default = 'Brownian')
%   Verbose: Verbosity level specifying how many temporary outputs should
%            be displayed (e.g., Command Window updates).
%
% OUTPUTS:
%   FitParams: Fit parameters for the CDF fit.  These will vary based on
%              'FitMethod'. FitParams(ii, :) will contain the fit
%              parameters for the fit to MSDStruct(ii).
%              (numel(MSDStruct)xNParameters array)
%   FitParamsSE: Standard errors for the CDF fit parameters 'FitParams'.
%                These will vary based on 'Method'.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Define default parameters if needed.
if (~exist('FitMethod', 'var') || isempty(FitMethod))
    FitMethod = 'WeightedLS';
end
if (~exist('NComponents', 'var') || isempty(NComponents))
    NComponents = 1;
end
if (~exist('DiffusionModel', 'var') || isempty(DiffusionModel))
    DiffusionModel = 'Brownian';
end
if (~exist('Verbose', 'var') || isempty(Verbose))
    Verbose = 0;
end

% Determine how many fit parameters are needed in the desired model.
NFitParams = 2 * NComponents;

% Fit the CDF data.
if (Verbose > 1)
    fprintf(['fitCDFOfJumps(): fitting CDF of jumps with FitMethod = ', ...
        '''%s''...\n'], FitMethod)
end
NFits = numel(MSDStruct);
FitParams = NaN(NFits, NFitParams);
FitParamsSE = NaN(NFits, NFitParams);
for ii = 1:NFits
    % Isolate some arrays needed for the fit.
    CDFOfJumps = double(MSDStruct(ii).CDFOfJumps);
    if isempty(CDFOfJumps)
        % This trajectory was too short so we can't make the fit.
        continue
    end
    FrameLagsAll = double(MSDStruct(ii).FrameLagsAll);
    NPoints = double(MSDStruct(ii).NPoints);
    LocVarianceSum = double(MSDStruct(ii).LocVarianceSum);
    SortedSquaredDisp = double(MSDStruct(ii).SortedSquaredDisp);
    
    % Fit the CDF of the jumps to the desired diffusion model.
    % NOTE: We're only allowing Brownian motion for now.
    switch lower(DiffusionModel)
        case 'brownian'
            if (nargout > 1)
                [ParamsHat, ParamsHatSE] = ...
                    smi_stat.DiffusionEstimator.fitCDFOfJumpsBrownian(...
                    SortedSquaredDisp, CDFOfJumps, ...
                    FrameLagsAll, NPoints, LocVarianceSum, NComponents, ...
                    [], FitMethod);
                FitParamsSE(ii, :) = ParamsHatSE.';
            else
                ParamsHat = ...
                    smi_stat.DiffusionEstimator.fitCDFOfJumpsBrownian(...
                    SortedSquaredDisp, CDFOfJumps, ...
                    FrameLagsAll, NPoints, LocVarianceSum, NComponents, ...
                    [], FitMethod);
            end
            FitParams(ii, :) = ParamsHat.';
        otherwise
            error('Unknown ''DiffusionModel'' = %s', DiffusionModel)
    end
end


end