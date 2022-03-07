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
%    % Data from some external source.
%    [x1, y1] = textread('9021_5.txt',  '%*u %u %u %*u', 'headerlines', 1);
%    SMD1.X = x1;
%    SMD1.Y = y1;
%    XY1 = [x1, y1];
%    [x2, y2] = textread('9021_10.txt', '%*u %u %u %*u', 'headerlines', 1);
%    SMD2.X = x2;
%    SMD2.Y = y2;
%    XY2 = [x2, y2];
%
%    PC.ROI = [0, 7400, 0, 6000];   % [x_min, x_max, y_min, y_max]
%    % These two numbers below are often the same, but the user can increase
%    % HistBinSize to make bigger internal pixels (or bigger internal histogram
%    % image bins) and so produce more smoothing, or decrease this quantity and
%    % attempt greater detail.
%    PC.PixelSize = 2.7559;         % Actual camera pixel size (nm).
%    PC.HistBinSize = 2.7559;       % Internal image pixel size (nm).
%
%    % Make a RoI structure (see also smi_cluster.ROITools).  This is typically
%    % used for invoking pair_correlation_ROIcombined with multiple ROIs, which
%    % are combined.  The size of each ROI need not be the same, although
%    % cleaner results will be produced if they are all the same size.
%    n_ROIs = 1;
%    RoI{1}.ROI = PC.ROI;
%    RoI{1}.X   = {XY1(:, 1), XY2(:, 1)};
%    RoI{1}.Y   = {XY1(:, 2), XY2(:, 2)};
%
%    % Typically, pair_correlation is used for comparing two images/ROIs, while
%    % pair_correlation_ROIcombined is used when combining the pair correlation
%    % results for several pairs of images/ROIs.  pair_correlation_Veatch is
%    % basically the original Sarah Veatch code left for comparison purposes.
%    % The results of all three routines in these examples should produce
%    % similar results (identical for the first two when n_ROIs = 1).
%    % NOTE: SMD structures or XY matrices can be provided as input.
%    PC.BaseName = '9021';
%    results_pcc  = PC.pair_correlation(SMD1, SMD2)
%    results_pcc  = PC.pair_correlation(XY1, XY2)
%    results_Rpcc = PC.pair_correlation_ROIcombined(2, n_ROIs, RoI)
%    results_Vpcc = PC.pair_correlation_Veatch(SMD1, SMD2, 'cross')
%    results_Vpcc = PC.pair_correlation_Veatch(XY1, XY2, 'cross')
%
%    PC.BaseName = '9021_5';
%    results_pac1  = PC.pair_correlation(SMD1)
%    results_pac1  = PC.pair_correlation(XY1)
%    results_Rpacc = PC.pair_correlation_ROIcombined(1, n_ROIs, RoI, 1)
%    results_Vpac1 = PC.pair_correlation_Veatch(SMD1,  [], 'auto')
%    results_Vpac1 = PC.pair_correlation_Veatch(XY1,  [], 'auto')
%
%    PC.BaseName = '9021_10';
%    results_pac2  = PC.pair_correlation(SMD2)
%    results_pac2  = PC.pair_correlation(XY2)
%    results_Rpac2 = PC.pair_correlation_ROIcombined(1, n_ROIs, RoI, 2)
%    results_Vpac2 = PC.pair_correlation_Veatch(SMD2, [], 'auto')
%    results_Vpac2 = PC.pair_correlation_Veatch(XY2, [], 'auto')

% =============================================================================
properties
% =============================================================================

   % If ROI is provided, it will be used, otherwise the xy_size will be
   % calculated using the (x_min, x_max, y_min, y_max) computed from the data
   % as
   %    xy_size = min(x_max - x_min, y_max - y_min)
   ROI = [];   % [x_min, x_max, x_max, y_max]   % nm

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
