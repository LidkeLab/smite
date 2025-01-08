function [M1, M2] = mandersSplitCoefs(SMD1, SMD2, SRZoom)
%mandersSplitCoefs computes the Mander's split coefficients between two SMDs.
% Manders' coefficients M1 and M2 measure colocalization by computing the
% intensity fraction of pixels from one channel that overlap with pixels from
% another.
% M1: intensity fraction of pixels in ch 1 that overlap with pixels in ch 2.
% M2: intensity fraction of pixels in ch 2 that overlap with pixels in ch 1.
%
% INPUTS:
%    SMD1, SMD2   single molecule data structures containing 2D localization
%                 coordinates in fields X, Y.
%    SRZoom       magnification factor for the SMD coordinates in order to get
%                 a better estimate of the Manders' Split Coefficients
%                 (default: 29)
%
% OUTPUTS:
%    M1, M2       M1 and M2 Manders' split coefficients.  Naming convention:
%                 https://imagej.net/imaging/colocalization-analysis
%
% CITATION:
%    E. M. M. Manders, F. J. Verbeek, and J. A. Aten
%    Measurement of co-localization of objects in dual-color confocal images.
%    Journal of Microscopy, Volume 169, Pt 3, March 1993, 375--382
%    https://imagej.net/media/manders.pdf

% Created by
%    Michael Wester (2024)

   if ~exist('SRZoom', 'var')
      SRZoom = 20;
   end

   % Make grayscale Gaussian blob images from the coordinates in SMD1 and SMD2.
   % 0s below are to omit the scalebar in the generated images.
   G1 = smi_vis.GenerateImages.grayscaleImage(SMD1, SRZoom, 0);
   G2 = smi_vis.GenerateImages.grayscaleImage(SMD2, SRZoom, 0);

   % Sum the intensities in channel 1 (N1) corresponding to colocalized pixels
   % in channel 2 that have intensities > 0, and similarly for channel 2.
   indx1 = G1(:) > 0;
   indx2 = G2(:) > 0;
   N1 = sum(G1(indx2));
   N2 = sum(G2(indx1));

   % Compute the total intensities in the two channels.
   Itotal1 = sum(G1(:));
   Itotal2 = sum(G2(:));

   % Manders' coefficients.
   M1 = N1 / Itotal1;
   M2 = N2 / Itotal2;

end
