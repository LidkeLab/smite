function success = unitTest()
% Test Fourier Ring Correlation (FRC) interface functions and provide examples
% of usage.
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
%    Install the Image Resolution software at the same level as sma-core-alpha.
%    This software is located at the URL above (see CITATION).  In startup.m,
%    add a path to
%       .../FRCresolution_software/matlabdistribution/FRCresolutionfunctions
%    where often ... = /Documents/MATLAB.  In the FRCresolutionfunctions, copy
%    smooth.m from the MATLAB Curve Fitting Toolbox into cfsmooth.m .  For
%    example, look in
%       MATLAB_DISTRIBUTION/toolbox/curvefit/curvefit/smooth.m
%    This is needed because DIPimage also has a smooth function which will
%    typically shadow MATLAB's smooth.

   SaveDir = smi_helpers.mkSMITETmpDir('unitTest', 'FRC');

   PixelSize = 100;   % nm
   SRZoom = 10;

   % Generate localizations and put them into an SMD structure.
   SIM = smi_sim.SimSMLM();
   SIM.SZ = 64;
   SIM.Rho = 50;
   SIM.NDatasets = 1;
   SIM.NFrames = 1000;
   SIM.simStar(16);
   SMD = SIM.SMD_Model;

   n_particles = numel(SMD.X);
   SMD.X_SE = (PixelSize/1000)*ones(n_particles, 1);
   SMD.Y_SE = (PixelSize/1000)*ones(n_particles, 1);

   %% ---------- uncorrected resolution calculation

   FRCc = smi_core.FRC();   % FRC class
   FRCc.PixelSize = PixelSize;
   FRCc.SRImageZoom = SRZoom;

   [res, ~, resH, resL] = FRCc.posToResolution(SMD);

   fprintf('resolution = %2.1f +- %2.2f [px]\n', ...
           res / SRZoom, (resL / SRZoom - resH / SRZoom)/2);
   fprintf('resolution = %2.1f +- %2.2f [nm]\n', ...
           res * PixelSize/SRZoom, (resL - resH)/2 * PixelSize/SRZoom);

   %% ---------- compute uncorrected FRC curve

   frc_curve = FRCc.posToFRC(SMD);
   [~, frc_curveA] = FRCc.posToResolution(SMD);

   figure();
   hold on
   qmax = 0.5 / (PixelSize/SRZoom);
   plot(linspace(0, qmax*sqrt(2), numel(frc_curve)), frc_curve,  'b-');
   plot(linspace(0, qmax*sqrt(2), numel(frc_curve)), frc_curveA, 'g-');
   xlim([0, qmax]);
   plot([0, qmax], [1/7, 1/7], 'r-');
   plot([0, qmax], [0, 0], 'k--');
   xlabel('spatial frequency (nm^{-1})');
   ylabel('FRC');
   title('Fourier Ring Correlation curve');
   legend({'posToFRC', 'posToResolution'}, 'Location', 'NorthEast');
   hold off
   saveas(gcf, fullfile(SaveDir, 'FRCcurve.png'));

   %% ---------- compute uncorrected and corrected FRC curves and resolutions
   % Correction removes spurious correlations and is the recommended method to
   % use for typical applications.

   try
   [res_corr, res_uncorr, Q, frc_curve_corr, frc_curve_uncorr] = ...
      FRCc.qCorrectionLocs(SMD);

   res1 = frctoresolution(frc_curve_uncorr, SMD.YSize*SRZoom);
   res2 = frctoresolution(frc_curve_corr,   SMD.YSize*SRZoom);
   fprintf('resolution = %2.1f, corrected = %2.1f [px], Q = %f\n', ...
           res_uncorr / SRZoom, res_corr / SRZoom, Q);
   fprintf('resolution = %2.1f, corrected = %2.1f [nm], Q = %f\n', ...
           res_uncorr * PixelSize/SRZoom, res_corr * PixelSize/SRZoom, Q);
   fprintf('frctoresolution = %2.1f, corrected = %2.1f [px]\n', res1, res2);

   figure();
   hold on
   qmax = 0.5 / (PixelSize/SRZoom);
   plot(linspace(0, qmax*sqrt(2), numel(frc_curve_corr)),   ...
        frc_curve_corr, 'b-');
   plot(linspace(0, qmax*sqrt(2), numel(frc_curve_uncorr)), ...
        frc_curve_uncorr, 'g-');
   xlim([0, qmax]);
   plot([0, qmax], [1/7, 1/7], 'r-');
   plot([0, qmax], [0, 0], 'k--');
   xlabel('spatial frequency (nm^{-1})');
   ylabel('FRC');
   title('Corrected Fourier Ring Correlation curve');
   legend({'corrected', 'uncorrected'}, 'Location', 'NorthEast');
   hold off
   saveas(gcf, fullfile(SaveDir, 'correctedFRCcurve.png'));
   end

   success = 1;

end
