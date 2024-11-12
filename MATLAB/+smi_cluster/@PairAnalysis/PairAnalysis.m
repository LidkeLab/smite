classdef PairAnalysis < handle

% PairAnalysis class written by Michael Wester
%    (10/13/2022) <wester@math.unm.edu>
% The New Mexico Center for the Spatiotemporal Modeling of Cell Signaling
% University of New Mexico Health Sciences Center
% Albuquerque, New Mexico, USA   87131
% Copyright (c) 2015-2023 by Michael J. Wester and Keith A. Lidke

% A collection of functions that interface with the main PairCorr routines.

% =============================================================================
properties
% =============================================================================

% =============================================================================
end % properties
% =============================================================================

% =============================================================================
methods
% =============================================================================

% =============================================================================
end % methods
% =============================================================================

% =============================================================================
methods(Static)
% =============================================================================

   n_ROIs_ALL = defineROIs2(Files1, Files2, Pixel2nm, Color, ROI_sizes, ...
                            ResultsDir, RegistrationNeeded, RegistrationViaDC);
   [results_pcc, resultsRC_pcc, results_c, results_cs, results_ls,   ...
    results_o1, results_o2] =                                        ...
      doAnalysis(n_ROIs, RoI, ROI_sizes, desc, particles, results_dir,      ...
                 options, PixelSize, HistBinSize, RmaxAxisLimit, E, minPts, ...
                 PlotNonOverlap, Color)
   results_bi = doBiStats(n_ROIs, RoI, desc, particles, results_dir, combined)
   results_c = doClustering(n_ROIs, RoI, desc, results_dir, plotting, ...
                            PixelSize, E, minPts)
   results_cs = ...
      doClusterSep2(n_ROIs, results_c, desc, particles, results_dir, plotting)
   results_ls = doLocSep2(n_ROIs, RoI, desc, particles, results_dir, plotting)
   [results_pcc, resultsRC_pcc] = ...
      doPairCorr(n_ROIs, RoI, ROI_sizes, desc, results_dir, combined, ...
                 plotting, HistBinSize, RmaxAxisLimit)
   results_o = doOverlap(n_ROIs, RoI, results_c, l12, desc, particles, ...
                         results_dir, PlotNonOverlap, Color, plotting)
   overlayBaGoLROIs(pathnameB, filesB, MAPNfile, ROI_sizes, SRImageZoom);
   doPlot2(n_ROIs, RoI, desc, particles, results_dir, Color, plotting)

% =============================================================================
end % methods(Static)
% =============================================================================
end % classdef
