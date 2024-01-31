classdef ClusterInterface < handle

% ClusterInterface class written by Michael Wester
%    (9/29/2022) <wester@math.unm.edu>
% The New Mexico Center for the Spatiotemporal Modeling of Cell Signaling
% University of New Mexico Health Sciences Center
% Albuquerque, New Mexico, USA   87131
% Copyright (c) 2015-2023 by Michael J. Wester and Keith A. Lidke

% A collection of functions that interface with the main Clustering routines.

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

   % Called by simpleROIcluster.
   combineResults(Files, analysis_dir, out_file)
   combinedStatistics1(SC, pathname, files, base_name, A_ROI, doHopkins)
   combinedStatistics2(SC, colors, line_type, pathname, files, ...
                       base_name, A_ROI, doHopkins)
   combineBaGoLROIs(pathnameR, filesR, pathnameB, filesB, MAPNfile, ...
                    keep_numbering)
   defineBaGoLROIs(pathnameR, filesR, pathnameB, filesB, MAPNfile, ...
                   OriginLLvsUL)
   defineROIs(pathname, files, Pixel2nm, RT, oneROI)
   [n_ROIs, RoI] = filterROIs(pathname, files, filter)
   % Called by singleConditionDriver.
   singleCondition(pathname, files, algorithm_range, E_range, ...
                   minPts_range, Pixel2nm, base_name, A_ROI,  ...
                   doHopkins, doSigmaActual, Alpha)

   % Called by plotROIDriver.
   plotROI(opt, pathnameC, filesC, pathnameB, filesB, PixelSize, SaveDir)
   plotROI1(opt, pathnameC, filesC, pathnameB, filesB, PixelSize, SaveDir)
   [SRIm] = genMAPNIm1(obj, ImFlag)
   OverlayImageCircle = genSRMAPNOverlay1(SMD, MAPN, XSize, YSize, ...
      PixelSize, SaveDir, Xstart, Ystart, RadiusScale, ScaleBarLength)

% =============================================================================
end % methods(Static)
% =============================================================================
end % classdef
