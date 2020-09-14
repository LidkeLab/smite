function [res_corr, res_uncorr, Q, frc_curve_corr, frc_curve_uncorr] = ...
         qCorrectionLocs(obj, SMR)
% qCorrectionLocs calculates Q-corrected FRC curves and resolutions.
%   
% INPUTS:
%    SMR             SMR structure with the following fields:
%       YSize        size (pixels) of the image
%       X,Y          (x, y) coordinates (pixels) of each localization [Nx1]
%       X_SE,Y_SE    standard errors of the estimation of the above   [NX1]
%       DatasetNum   dataset number of each localization              [Nx1]
%       FrameNum     frame number of each localization                [Nx1]
%       NFrames      number of frames in each dataset
%    obj:
%       PixelSize       nm per pixel
%       SRImageZoom     image magnification factor
%       Repeats         number of times the FRC curve is computed for averaging
%
% OUTPUTS:
%    res_corr         resolution after correction for spurious correlations
%    res_uncorr       resolution value without correction
%    Q                estimate for the number of times an emitter is localized
%                        on average assuming Poisson statistics for the
%                        localizations per emitter
%    frc_curve_corr   FRC curve corrected for spurious correlations
%    frc_curve_uncorr FRC curve without correction
% Note: frc_curve is the FRC curve vs. spatial frequency (1/nm)
%
% REQUIRES:
%    DIPlib Image Resolution add-on
%    Curve Fitting Toolbox
%    Parallel Processing Toolbox
%    NVidia GPU
%  
% CITATION:
%    http://www.diplib.org/add-ons
%    Image Resolution, Reference: R.P.J. Nieuwenhuizen, K.A. Lidke, M. Bates,
%    D. Leyton Puig, D. Grünwald, S. Stallinga, B. Rieger, Measuring Image
%    Resolution in Optical Nanoscopy, Nature Methods, 10(6):557-562, 2013.
%    Also, see FRC_unitTest.m .

   % Compute absolute frame numbers as if there is just one dataset.
   FrameAbs = double(SMR.NFrames)*(SMR.DatasetNum - 1) + double(SMR.FrameNum);
   SRImageZoom = obj.SRImageZoom;
   % Number of time blocks into which the dataset is split.  2 blocks is the
   % minimum accepted amount.
   blocks = 50;
   reps = obj.Repeats;
   sig = sqrt(SMR.X_SE.^2 + SMR.Y_SE.^2);
   meansig = mean(sig);
   stdsig  = std(sig);
   SR_pixelsize = obj.PixelSize; % nm
   floorcor = true;
   show_frc = false;

   [res_corr, res_uncorr, Q, frc_curve_corr, frc_curve_uncorr] = ...
      qcorrection_locs([SMR.X, SMR.Y, FrameAbs], SMR.YSize*SRImageZoom, ...
                        SRImageZoom, blocks, reps, meansig, stdsig,     ...
                        SR_pixelsize, floorcor, show_frc);

end
