function [MLEParams, MLEParamsSE] = ...
    mleOfJumps(MSDStruct, DiffusionModel, Verbose)
%mleOfJumps finds the MLE from the likelihood of observed jumps.
% This method will make a maximum likelihood estimate of the diffusion
% constant(s) (and possibly population ratios) based on a model for the
% likelihood of observed jumps.
%
% INPUTS:
%   MSDStruct: A structure array of MSD data as output from computeMSD()
%              (see computeMSD() for more details).
%   DiffusionModel: A string specifying the diffusion model to fit to the
%                   MSD. See options in DiffusionEstimator class property
%                   'DiffusionModel'. (Default = 'Brownian2C')
%   Verbose: Verbosity level specifying how many temporary outputs should
%            be displayed (e.g., Command Window updates). (Default = 0)
%
% OUTPUTS:
%   MLEParams: MLE for the desired model parameters..  These will vary 
%              based on 'FitMethod'. MLEParams(ii, :) will contain the MLEs 
%              corresponding to the data in MSDStruct(ii).
%              (numel(MSDStruct)xNParameters array)
%   MLEParamsSE: Standard errors for the MLEs in 'MLEParams'. These will
%                vary based on 'Method'.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Define default parameters if needed.
if (~exist('DiffusionModel', 'var') || isempty(DiffusionModel))
    DiffusionModel = 'Brownian2C';
end
if (~exist('Verbose', 'var') || isempty(Verbose))
    Verbose = 0;
end

% Determine how many components and fit parameters are used in the desired
% model.
switch lower(DiffusionModel)
    case 'brownian1c'
        NComponents = 1;
    case 'brownian2c'
        NComponents = 2;
    otherwise
        error('Unknown ''DiffusionModel'' = %s', DiffusionModel)
end
NFitParams = 2*NComponents - 1;

% Compute the MLE.
if (Verbose > 1)
    fprintf('mleOfJumps(): computing MLE from observed jumps...\n')
end
NFits = numel(MSDStruct);
MLEParams = NaN(NFits, NFitParams);
MLEParamsSE = NaN(NFits, NFitParams);
for ii = 1:NFits
    % Isolate some arrays needed for the fit.
    SquaredDisplacement = double(MSDStruct(ii).SquaredDisplacement);
    if isempty(SquaredDisplacement)
        % This trajectory was too short and no displacements were measured.
        continue
    end
    FrameLagsAll = double(MSDStruct(ii).FrameLagsAll);
    LocVarianceSum = double(MSDStruct(ii).LocVarianceSum);
    
    % Compute the MLE based on the desired model.
    % NOTE: We're only allowing Brownian motion for now.
    if (nargout > 1)
        [ParamsHat, ParamsHatSE] = ...
            smi_stat.DiffusionEstimator.mleOfJumpsBrownian(...
            SquaredDisplacement, FrameLagsAll, LocVarianceSum, NComponents);
        MLEParamsSE(ii, :) = ParamsHatSE.';
    else
        ParamsHat = smi_stat.DiffusionEstimator.mleOfJumpsBrownian(...
            SquaredDisplacement, FrameLagsAll, LocVarianceSum, NComponents);
    end
    MLEParams(ii, :) = ParamsHat.';
end


end