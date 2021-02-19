function [FitParams, FitParamsSE] = fitMSD(MSDStruct, FitMethod, Verbose)
%fitMSD fits mean squared displacement data.
% This method will fit a mean squared displacement plot by the method
% specified by 'Method'.
%
% INPUTS:
%   MSDStruct: A structure array of MSD data as output from computeMSD()
%              (see computeMSD() for more details).
%   FitMethod: A string specifying the fit method. (Default = 'weightedLS')
%   Verbose: Verbosity level specifying how many temporary outputs should
%            be displayed (e.g., Command Window updates).
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
if (~exist('Method', 'var') || isempty(FitMethod))
    FitMethod = 'weightedLS';
end
if (~exist('Verbose', 'var') || isempty(Verbose))
    Verbose = 0;
end

% Fit the MSD data.
if (Verbose > 1)
    fprintf('fitMSD(): fitting MSD with FitMethod = ''%s''...\n', ...
        FitMethod)
end
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
    switch FitMethod
        case 'weightedLS'
            % Fit the MSD using weighted least squares.
            % NOTE: The weight array is (as I've defined it here) 
            %       proportional to the reciprocal of the CRLB of each 
            %       point in an MSD plot. The proportionality constant is 
            %       dropped because it doesn't affect the weighted fit.
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