classdef PairCorrelation < handle

% PairCorrelation class written by Michael Wester, Keith Lidke,
%    and others as noted internally (2/26/2018) <wester@math.unm.edu>
% The New Mexico Center for the Spatiotemporal Modeling of Cell Signaling
% University of New Mexico Health Sciences Center
% Albuquerque, New Mexico, USA   87131
% Copyright (c) 2015-2021 by Michael J. Wester and Keith A. Lidke
%
% PairCorrelation performs pair correlation (auto- and cross-correlation) on
% user selected rectangular regions of interest (ROIs).  Input units are
% assumed to be in units of nm.  Functional dependencies:
%
%    pair_correlation             <- get_autocorr, get_crosscorr, pc_GaussFit
%    pair_correlation_Veatch      <- get_autocorr, get_crosscorr, pc_GaussFit
%    pair_correlation_ROIcombined <- get_corr
%
% pair_correlation (and pair_correlation_Veatch) operates on a single pair of
% images (pair_correlation_Veatch is basically the original code written by
% Sarah Veatch).  pair_correlation_ROIcombined combines pair correlation on a
% series of images/ROIs.
%
% Example main program:
%
%    PC = smi_cluster.PairCorrelation();
%    PC.ResultsDir = 'Results';
%
%    [x, y] = textread('9021_5.txt',  '%*u %u %u %*u', 'headerlines', 1);
%    XY_5 =  [x, y];
%    [x, y] = textread('9021_10.txt', '%*u %u %u %*u', 'headerlines', 1);
%    XY_10 = [x, y];
%    PC.ROI = [0, 7400, 0, 6000];   % [xmin, xmax, ymin, ymax]
%    PC.HistBinSize = 2.7559;
%
%    % See also smi_cluster.ROITools
%    n_ROIs = 1;
%    RoI{1}.ROI = PC.ROI;
%    RoI{1}.X   = {XY_5(:, 1), XY_10(:, 1)};
%    RoI{1}.Y   = {XY_5(:, 2), XY_10(:, 2)};
%    
%    PC.BaseName = '9021';
%    results_pcc  = PC.pair_correlation(XY_5, XY_10)
%    results_Vpcc = PC.pair_correlation_Veatch(XY_5, XY_10, 'cross')
%    results_Rpcc = PC.pair_correlation_ROIcombined(2, n_ROIs, RoI)
%
%    pc.BaseName = '9021_5';
%    results_pac1  = PC.pair_correlation(XY_5)
%    results_Vpac1 = PC.pair_correlation_Veatch(XY_5,  [], 'auto')
%    results_Rpcc  = PC.pair_correlation_ROIcombined(1, n_ROIs, RoI, 1)
%
%    PC.BaseName = '9021_10';
%    results_pac2  = PC.pair_correlation(XY_10)
%    results_Vpac2 = PC.pair_correlation_Veatch(XY_10, [], 'auto')
%    results_Rpcc  = PC.pair_correlation_ROIcombined(1, n_ROIs, RoI, 2)

% =============================================================================
properties
% =============================================================================

   % If ROI is provided, it will be used, otherwise the xy_size will be
   % calculated using the (x_min, x_max, y_min, y_max) computed from the data
   % as
   %    xy_size = min(x_max - x_min, y_max - y_min)
   ROI = [];   % [x_min, x_max, x_max, y_max]

   BaseName = ''; % descriptive name for the result files.
   Font_props = {'FontSize', 15, 'FontWeight', 'bold'};
   % If Fig_ext is empty, plot to the screen and save as .fig; otherwise, do
   % not plot to the screen, but save as .fig AND .Fig_ext .
   Fig_ext = 'png';
   Lines   = true;  % If true, plot lines rather than points for g(r) vs. r
   ResultsDir = '.';   % Directory to store results.

   % Pixel size (nm) used for the internal images to be correlated.
   HistBinSize = 100;
   PixelSize = 100;   % Actual camera pixel size (nm).
   Rmax_axis = -1;   % Sets plotting limit if > 0 (nm)
   % Default fit model for pair_correlation_Veatch.  This can also be supplied
   % as the optional last (4th) argument to this function.  Possible choices:
   %    'exponential_and_gaussian', 'exponential_and_cosine', 'exponential'
   Veatch_fit = 'exponential_and_gaussian';
   % Note that pair_correlation and pair_correlation_ROIcombined use a 2D
   % Gaussian fit model.

   Verbose = 1;   % verbosity level

% =============================================================================
end % properties

methods
% =============================================================================

% Constructor.  SMF is an optional argument.
function obj = PairCorrelation(SMF)

   if ~exist('SMF', 'var')
      SMF = smi_core.SingleMoleculeFitting();
   end
   obj.ResultsDir  = SMF.Data.ResultsDir;
   obj.HistBinSize = SMF.Data.PixelSize;
   obj.PixelSize   = SMF.Data.PixelSize;

end

% =============================================================================
end % methods

methods(Static)
% =============================================================================

   [G, r, g, dg, mask, rmax] = get_autocorr(I1, mask, rmax, flag)
   [C, r, c, dc, mask, rmax] = get_crosscorr(I1, I2, mask, rmax, flag)
   [C, r, c, dc, rmax] = get_corr(n_ROIs, rmax, II1, II2)
   [estimates, errors, model] = pc_GaussFit(r,g_r,rmax,rho)

% =============================================================================
end % methods(Static)
% =============================================================================
end % classdef
