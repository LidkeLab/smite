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
%              'FitMethod'. FitParams(ii, :) will contain the fit
%              parameters for the fit to MSDStruct(ii). 
%              (numel(MSDStruct)xNParameters array)
%   FitParamsSE: Standard errors for the MSD fit parameters 'FitParams'.
%                These will vary based on 'Method'.
%
% REQUIRES:
%   Curve Fitting Toolbox (for fit() and associated methods)

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Define default parameters if needed.
if (~exist('FitMethod', 'var') || isempty(FitMethod))
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
NFits = numel(MSDStruct);
FitParams = NaN(NFits, 2);
FitParamsSE = NaN(NFits, 2);
for ii = 1:NFits
    % If the MSD has enough points to be "useful", try to fit it.
    FrameLags = double(MSDStruct(ii).FrameLags);
    NFrames = numel(FrameLags);
    NPoints = double(MSDStruct(ii).NPoints);
    MSD = double(MSDStruct(ii).MSD);
    if (NFrames < 3)
        continue
    end
    switch FitMethod
        case 'LS'
            % Fit the MSD using linear least squares.
            [BetaHat, BetaHatSE] = ...
                smi_stat.leastSquaresFit(FrameLags, MSD);
            FitParams(ii, :) = BetaHat.';
            FitParamsSE(ii, :) = BetaHatSE.';
        case 'weightedLS'
            % Fit the MSD using weighted least squares.
            % NOTE: The weight array is (as I've defined it here) 
            %       proportional to the reciprocal of the CRLB of each 
            %       point in an MSD plot. The proportionality constant is 
            %       dropped because it doesn't affect the weighted fit.
            Weights = NPoints ./ (FrameLags.^2);
            [BetaHat, BetaHatSE] = ...
                smi_stat.leastSquaresFit(FrameLags, MSD, Weights);
            FitParams(ii, :) = BetaHat.';
            FitParamsSE(ii, :) = BetaHatSE.';
        otherwise
            error('Unknown ''FitMethod'' = %s', FitMethod)
    end
end


end