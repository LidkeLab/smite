classdef FRC < handle

% Fourier Ring Correlation (FRC) class for computing average image resolution.
% FRC is a measure of the average resolution over a super-resolution image.  It
% works by dividing the set of single-emitter localizations in the
% super-resolution image into two statistically independent subsets.  The
% Fourier transforms of subimages generated from each of these subsets are then
% statistically correlated over pixels on the perimeters of circles of constant
% spatial frequency.  The image resolution is defined as the inverse of the
% spatial frequency when the FRC curve drops below a threshold, taken to be
% 1/7 in the citation below.  Spurious correlations (for example, due to
% repeated photoactivation of the same emitter) are removed by estimating the
% number of times an emitter is localized on average (Q) assuming Poisson
% statistics.
%
% REQUIRES:
%    DIPlib Image Resolution add-on
%    Curve Fitting Toolbox (needed by qCorrectLocs)
%    Parallel Processing Toolbox
%    NVidia GPU
%
% CITATION:
%    http://www.diplib.org/add-ons
%    Image Resolution, Reference: R.P.J. Nieuwenhuizen, K.A. Lidke, M. Bates,
%    D. Leyton Puig, D. Grünwald, S. Stallinga, B. Rieger, Measuring Image
%    Resolution in Optical Nanoscopy, Nature Methods, 10(6):557-562, 2013.
% NOTE:
%    Install the Image Resolution software at the same level as smite.
%    This software is located at the URL above (see CITATION).  In startup.m,
%    add a path to
%       .../FRCresolution_software/matlabdistribution/FRCresolutionfunctions
%    where often ... = /Documents/MATLAB.  In the FRCresolutionfunctions, copy
%    smooth.m from the MATLAB Curve Fitting Toolbox into cfsmooth.m .  For
%    example, look in
%       MATLAB_DISTRIBUTION/toolbox/curvefit/curvefit/smooth.m
%    This is needed because DIPimage also has a smooth function which will
%    typically shadow MATLAB's smooth. 

% =============================================================================
properties

   PixelSize   = 100;   % nm per pixel
   SRImageZoom = 10;    % image magnification factor
   Repeats     = 1;     % number of times FRC curve is computed for averaging

end % properties
% =============================================================================

% =============================================================================
methods

   frc_curve = posToFRC(obj, SMR)
   [resolution, frc_curve, resolution_hi, resolution_lo, resolution_t] = ...
      posToResolution(obj, SMR)
   [res_corr, res_uncorr, Q, frc_curve_corr, frc_curve_uncorr] = ...
      qCorrectionLocs(obj, SMR)

   % Constructor.
   function obj = FRC(SMF)
   % SMF values, if provided, can override class properties.

      if exist('SMF', 'var')
         if ~isempty(FRC.Data.PixelSize)
            obj.PixelSize = 1000 * FRC.Data.PixelSize;   % convert um to nm
         end
      end

   end

end % methods
% =============================================================================

% =============================================================================
methods(Static)

   success = unitTest()

end % methods(Static)
% =============================================================================

end % classdef FRC
