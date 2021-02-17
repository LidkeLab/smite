function [FitParams, FitParamsSE] = fitMSD(MSDStruct, Method)
%fitMSD fits mean squared displacement data.
% This method will fit a mean squared displacement plot by the method
% specified by 'Method'.
%
% INPUTS:
%   MSDStruct: A structure array of MSD data as output from computeMSD()
%              (see computeMSD() for more details).
%   Method: A string specifying the fit method. (Default = 'weightedLS')
%
% OUTPUTS:
%   FitParams: Fit parameters for the MSD fit.  These will vary based on
%              'Method'.
%   FitParamsSE: Standard errors for the MSD fit parameters 'FitParams'.
%                These will vary based on 'Method'.
%
% REQUIRES:
%   Curve Fitting Toolbox

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Define default parameters if needed.
if (~exist('Method', 'var') || isempty(Method))
    Method = 'weightedLS';
end

% Fit the MSD data.
NTraj = numel(MSDStruct);
FitParams = NaN(NTraj, 2);
FitParamsSE = NaN(NTraj, 2);
for ii = 1:NTraj
    % If the MSD has enough points to be "useful", try to fit it.
    FrameLags = double(MSDStruct(ii).FrameLags);
    NPoints = double(MSDStruct(ii).NPoints);
    MSD = double(MSDStruct(ii).MSD);
    if (numel(FrameLags) < 3)
        continue
    end
    switch Method
        case 'weightedLS'
            % Fit the MSD using weighted least squares.
            WeightArray = NPoints ./ (FrameLags.^2);
            FitResults = fit(FrameLags, MSD, 'poly1', ...
                'Weights', WeightArray);
            
            % Extract the parameters of interest from the fit results.
            % NOTE: I could just specify the ~0.68 confidence level and
            %       remove the factor of 1.96, but these numbers are more
            %       round....
            FitParams(ii, :) = coeffvalues(FitResults);
            FitParamsCI = confint(FitResults, 0.95);
            FitParamsSE(ii, :) = ...
                (FitParamsCI(2, :) - FitParams(ii, :)) / 1.96;
    end
end


end