function [FitParams, FitParamsSE] = ...
    fitMSD(MSDStruct, DiffusionModel, FitMethod, Verbose)
%fitMSD fits mean squared displacement data.
% This method will fit a mean squared displacement plot by the method
% specified by 'Method'.
%
% INPUTS:
%   MSDStruct: A structure array of MSD data as output from computeMSD()
%              (see computeMSD() for more details).
%   DiffusionModel: A string specifying the diffusion model to fit to the
%                   MSD. See options in DiffusionEstimator class property
%                   'DiffusionModel'. (Default = 'Brownian')
%   FitMethod: A string specifying the fit method. (Default = 'WeightedLS')
%   Verbose: Verbosity level specifying how many temporary outputs should
%            be displayed (e.g., Command Window updates).
%
% OUTPUTS:
%   FitParams: Fit parameters for the MSD fit.  These will vary based on
%              'FitMethod'. FitParams(ii, :) will contain the fit
%              parameters for the fit to MSDStruct(ii). 
%              (numel(MSDStruct)xNParameters array)
%   FitParamsSE: Standard errors for the MSD fit parameters 'FitParams'.
%                These will vary based on 'Method'.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Define default parameters if needed.
if (~exist('DiffusionModel', 'var') || isempty(DiffusionModel))
    DiffusionModel = 'Brownian';
end
if (~exist('FitMethod', 'var') || isempty(FitMethod))
    FitMethod = 'WeightedLS';
end
if (~exist('Verbose', 'var') || isempty(Verbose))
    Verbose = 0;
end

% Fit the MSD data.
if (Verbose > 1)
    fprintf('fitMSD(): fitting MSD with FitMethod = ''%s''...\n', ...
        FitMethod)
end
NFits = numel(MSDStruct);
FitParams = NaN(NFits, 2);
FitParamsSE = NaN(NFits, 2);
for ii = 1:NFits
    % Make sure the MSD has enough points to make a useful fit.
    FrameLags = double(MSDStruct(ii).FrameLags);
    NFrames = numel(FrameLags);
    if (NFrames < 2)
        continue
    end
    
    % Fit the MSD to the desired diffusion model.
    NPoints = double(MSDStruct(ii).NPoints);
    MSD = double(MSDStruct(ii).MSD);
    switch DiffusionModel
        case {'Brownian', 'brownian'}
            [FitParams, FitParamsSE] = ...
                smi_stat.DiffusionEstimator.fitMSDBrownian(...
                FrameLags, MSD, NPoints, FitMethod);
        otherwise
            error('Unknown ''DiffusionModel'' = %s', DiffusionModel)
    end
end


end