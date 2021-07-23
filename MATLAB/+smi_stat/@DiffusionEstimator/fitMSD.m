function [FitParams, FitParamsSE] = ...
    fitMSD(MSDStruct, FitMethod, NFitPoints, DiffusionModel, Verbose)
%fitMSD fits mean squared displacement data.
% This method will fit a mean squared displacement plot by the method
% specified by 'FitMethod'.
%
% INPUTS:
%   MSDStruct: A structure array of MSD data as output from computeMSD()
%              (see computeMSD() for more details).
%   FitMethod: A string specifying the fit method. (Default = 'WeightedLS')
%   NFitPoints: Number of points in the MSD to be fit. (Default = 5)
%   DiffusionModel: A string specifying the diffusion model to fit to the
%                   MSD. See options in DiffusionEstimator class property
%                   'DiffusionModel'. (Default = 'Brownian')
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
if (~exist('FitMethod', 'var') || isempty(FitMethod))
    FitMethod = 'WeightedLS';
end
if (~exist('NFitPoints', 'var') || isempty(NFitPoints))
    NFitPoints = 5;
end
if (~exist('DiffusionModel', 'var') || isempty(DiffusionModel))
    DiffusionModel = 'Brownian';
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
    MaxFitPoints = numel(FrameLags);
    FitPointsIndices = 1:min(NFitPoints, MaxFitPoints);
    FrameLags = FrameLags(FitPointsIndices);
    if (numel(FrameLags) < 2)
        % This trajectory is too short so we won't bother fitting it.
        continue
    end
    
    % Fit the MSD to the desired diffusion model.
    NPoints = double(MSDStruct(ii).NPoints(FitPointsIndices));
    MSD = double(MSDStruct(ii).MSD(FitPointsIndices));
    switch lower(DiffusionModel)
        case 'brownian1c'
            [FitParams(ii, :), FitParamsSE(ii, :)] = ...
                smi_stat.DiffusionEstimator.fitMSDBrownian(...
                FrameLags, MSD, NPoints, FitMethod);
        case 'brownian2c'
            error(['The two component model ''Brownian2C'' can only ', ...
                'be used when obj.FitTarget = ''CDFOfJumps'''])
        otherwise
            error('Unknown ''DiffusionModel'' = %s', DiffusionModel)
    end
end


end