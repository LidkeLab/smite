function plotROI(opt, pathnameC, filesC, pathnameB, filesB, PixelSize, SaveDir)
% Plot dot, Gaussian or circle images of SMD/MAPN coordinates per ROI per cell.
%
% INPUTS:
%    opt         data characteristics and types of plots to produce if true
%       SR          SR Results file
%       BaGoL       BaGoL Results file (BGL.SMD)
%       MAPN        BaGoL MAPN file
%       Dot         Dot plot
%       Gaussian    Gaussian plot
%       Circle      Circle plot (BaGoL Results file: BGL.SMD + BGL.MAPN)
%       Boundary    Include ROI boundaries
%       Cluster     Include ROI clusters
%       NoSave      Do not save outputs
%    pathnameC   path to filesC
%    filesC      cluster data per ROI per cell for a single condition; this
%                will be a single *_results.mat file
%    pathnameB   path to filesB
%    filesB      Results files (SR, BaGoL or MAPN) defining SMD-like structures,
%                typically representing several cell images collected under a
%                single experimental condition and clustered together
%    PixelSize   conversion factor from pixels to nm
%    SaveDir     directory to which the plots produced are saved
%
% OUTPUTS:
%    The various plots opted for in the input will be produced for each ROI in
%    each fileB.

% Created by
%    Michael J. Wester (2022)

   % If zero, ROIs in simpleROIcluster were chosen using plotted points, but
   % if > 0, ROIs were chosen using GaussianIm = true so following the DIPimage
   % convention---the value should then be YSize * PixelSize.
   % Changes (9/20/21) in simpleROIcluster seem to have made this fix obsolete,
   % but will retain the code for new just in case they are needed once again.
   % Certainly useful when plotting boundaries/clusters on top of images.
   GaussianImageKludge = 256*PixelSize;
   % Boundary shrink factor (0 = convex hull, 1 = as concave as possible,
   % 0.5 = MATLAB default)
   ShrinkFactor = 0.5;
   % Zoom factor for SR images
   Zoom = 20;
   ScaleBarLength = 500; % nm
   ScaleBarWidth  = 100; % nm

   CI = smi_cluster.ClusterInterface();

   dataC = load(fullfile(pathnameC, filesC{1}));
   for i = 1 : numel(filesB)
      fileB = filesB{i};
      dataB = load(fullfile(pathnameB, fileB));
      if opt.SR
         SMD = dataB.SMD;
      end
      if opt.MAPN
         MAPN = dataB.MAPN;
      end
      % Remove extraneous material from the filenames.
      short = regexprep(fileB, '.mat$', '');
      short = regexprep(short, '_ResultsStruct$', '');
      short = regexprep(short, '_Results$', '');
      short = regexprep(short, '^BaGoL_Results_', '');
      shrt  = regexprep(short, '_', '\\_');

      % Display a quick summary of the upcoming analysis.
      if opt.SR
         fprintf('%s: SMD = %d\n', short, numel(SMD.X));
      elseif opt.BaGoL || opt.Circle
         BGL = dataB.BGL;
         fprintf('%s: SMD = %d, MAPN = %d\n', ...
                 short, numel(BGL.SMD.X), numel(BGL.MAPN.X));
      elseif opt.MAPN
         fprintf('%s: MAPN = %d\n', short, numel(MAPN.X));
      end

      % SR units conversion.
      if opt.SR
         SMD.X = SMD.X * PixelSize;
         SMD.Y = GaussianImageKludge - SMD.Y * PixelSize;
         SMD.X_SE = SMD.X_SE * PixelSize;
         SMD.Y_SE = SMD.Y_SE * PixelSize;
      end

      % Plot ROIs individually.
      if opt.SR
         SMDsave = SMD;
      elseif opt.BaGoL
         SMDsave = BGL.SMD;
      elseif opt.MAPN
         %SMDsave = BGL.MAPN;
         SMDsave = MAPN;
      elseif opt.Circle
         SMDsave  = BGL.SMD;
         MAPNsave = BGL.MAPN;
      end
      j_ROI = sum(dataC.n_ROIs(1 : i - 1));
      for j = 1 : dataC.n_ROIs(i)
         j_ROI = j_ROI + 1;
         ROI = dataC.RoI{i}{j}.ROI;
         SMD = SMDsave;

         if opt.Boundary
            xExtra = 0.1 * (ROI(2) - ROI(1));
            yExtra = 0.1 * (ROI(4) - ROI(3));
         else
            xExtra = 0;
            yExtra = 0;
         end

         if opt.Dot
            plot(SMD.X, SMD.Y, 'k.');
            hold on
         elseif opt.Gaussian
            if opt.MAPN
               BGL.MAPN = SMD;
            end
            % Gaussian image plot of SR localizations.
            BGL.PixelSize = 1;
            %if GaussianImageKludge > 0
            %   BGL.MAPN.Y = GaussianImageKludge - BGL.MAPN.Y;
            %end
if ~opt.SR
            GIK = GaussianImageKludge;
            indx = ROI(1) - xExtra <= SMD.X & SMD.X <= ROI(2) + yExtra & ...
                   ROI(3) - yExtra <= GIK - SMD.Y &                      ...
                   GIK - SMD.Y <= ROI(4) + yExtra;
            SMD.X = SMD.X(indx);
            SMD.Y = SMD.Y(indx);
            SMD.X_SE = SMD.X_SE(indx);
            SMD.Y_SE = SMD.Y_SE(indx);
end
            if opt.BaGoL
               BGL.SMD = SMD;
               MapIm = CI.genMAPNIm1(BGL, 2);
            elseif opt.MAPN
               BGL.PImageSize = dataB.PImageSize;
               BGL.XStart = dataB.XStart;
               BGL.YStart = dataB.YStart;
               if GaussianImageKludge > 0
                  BGL.MAPN.Y = GaussianImageKludge - BGL.MAPN.Y;
               end
               MapIm = CI.genMAPNIm1(BGL, 1);
            elseif opt.SR && j == 1
               %MapIm = BGL.makeIm(SMD, 256*PixelSize, PixelSize);
               SMD_SR = SMD;
               SMD_SR.X = SMD_SR.X ./ PixelSize;
               SMD_SR.Y = SMD_SR.Y ./ PixelSize;
               SMD_SR.X_SE = SMD_SR.X_SE ./ PixelSize;
               SMD_SR.Y_SE = SMD_SR.Y_SE ./ PixelSize;
               MapIm = smi_vis.GenerateImages.gaussianImage(SMD_SR, Zoom, ...
                  ScaleBarLength/1000*Zoom);
               imshow(MapIm, hot);
               hold on
            end
            if ~opt.SR
               MapIm = smi.BaGoL.scaleIm(MapIm, 98);
            end
         end

         % Gaussian image plot of SR localizations.
         ScaleBarLength = 100;  % nm
         ScaleBarWidth  = 25;   % nm
         if opt.Gaussian
            if j == 1
               fprintf('ScaleBar = %d nm\n', ScaleBarLength);
               fprintf('ROI (%d):', dataC.n_ROIs(i));
            end
            fprintf(' %d', j);
            if ~opt.SR
               % Add in a scale bar (lower right).
               Xoffset = 100;
               Yoffset = 100;
               Xstart = ROI(2) - Xoffset - ScaleBarLength;
               Ystart = ...
                  GaussianImageKludge - (ROI(3) + Yoffset + ScaleBarWidth);
               X = round(Xstart) : round(Xstart + ScaleBarLength);
               Y = round(Ystart) : round(Ystart + ScaleBarWidth);
               MapIm(Y, X) = 255;
               imshow(MapIm, hot);
            end

            if opt.SR && j > 1
               imshow(MapIm, hot);
            end
            hold on
         elseif opt.Circle
            BGL.PixelSize = 1;
            SMD = SMDsave;
            GIK = GaussianImageKludge;
            indx = ROI(1) <= SMD.X & SMD.X <= ROI(2) & ...
                   ROI(3) <= GIK - SMD.Y & GIK - SMD.Y <= ROI(4);
            SMD.X = SMD.X(indx) - ROI(1);
            SMD.Y = SMD.Y(indx) - (GIK - ROI(4));
            SMD.X_SE = SMD.X_SE(indx);
            SMD.Y_SE = SMD.Y_SE(indx);
            BGL.SMD = SMD;

            MAPN = MAPNsave;
            GIK = GaussianImageKludge;
            indx = ROI(1) <= MAPN.X & MAPN.X <= ROI(2) & ...
                   ROI(3) <= GIK - MAPN.Y & GIK - MAPN.Y <= ROI(4);
            MAPN.X = MAPN.X(indx) - ROI(1);
            MAPN.Y = MAPN.Y(indx) - (GIK - ROI(4));
            MAPN.X_SE = MAPN.X_SE(indx);
            MAPN.Y_SE = MAPN.Y_SE(indx);
            BGL.MAPN = MAPN;

            ROISize = ROI(2) - ROI(1);
            MapIm = CI.genSRMAPNOverlay1(BGL.SMD, BGL.MAPN, ROISize,      ...
                                         ROISize, 1, SaveDir, BGL.XStart, ...
                                         BGL.YStart, 1, ScaleBarLength);

            % Gaussian image plot of SR (green) + MAPN (magenta) localizations.
            imshow(MapIm, hot);
            hold on
         end

         % Plot ROI boundaries.
         if opt.Boundary
            if opt.Circle
               y = [ROI(1), ROI(2), ROI(2), ROI(1), ROI(1)];
               x = [ROI(3), ROI(3), ROI(4), ROI(4), ROI(3)];
               %if GaussianImageKludge > 0
               %   y = GaussianImageKludge - y;
               %end
            else
               x = [ROI(1), ROI(2), ROI(2), ROI(1), ROI(1)];
               y = [ROI(3), ROI(3), ROI(4), ROI(4), ROI(3)];
               if GaussianImageKludge > 0
                  y = GaussianImageKludge - y;
               end
            end
            if opt.SR && opt.Gaussian
               x = x ./ PixelSize * Zoom;
               y = y ./ PixelSize * Zoom;
            end
            plot(x, y, 'g-', 'LineWidth', 2);
         end

         if  yExtra > 0
            % Add a label designating the ROI number.
            x_label = (ROI(1) + ROI(2))/2;
            y_label = ROI(4) + yExtra/2;
            if GaussianImageKludge > 0
               y_label = GaussianImageKludge - y_label;
            end
            if opt.SR && opt.Gaussian
               x_label = x_label ./ PixelSize * Zoom;
               y_label = y_label ./ PixelSize * Zoom;
            end
            text(x_label, y_label, sprintf('%d', j), 'Color', 'green');
         end

         if opt.Cluster
            % Plot cluster boundaries.
            results = dataC.results{j_ROI};
            for l = 1 : results.nC
               x = results.XY(:, 1);
               y = results.XY(:, 2);
               if GaussianImageKludge > 0
                  y = GaussianImageKludge - y;
               end
               if opt.SR && opt.Gaussian
                  x = x ./ PixelSize * Zoom;
                  y = y ./ PixelSize * Zoom;
               end
               C = results.C;
               if numel(C{l}) <= 2
                  plot(x(C{l}), y(C{l}), 'm.-', 'LineWidth', 2, ...
                       'MarkerSize', 12);
               else
                  % Determine cluster boundary indices from the cluster
                  % contents.
                  k = boundary(double(x(C{l})), double(y(C{l})), ShrinkFactor);
                  k = C{l}(k);
                  k = [k, k(1)];
                  plot(x(k), y(k), 'm.-', 'LineWidth', 2, 'MarkerSize', 12);
               end
            end % for l (cluster #)
         end % if opt.Cluster

         % Circle plots already are the full ROI.
         if ~opt.Circle
            if opt.SR
               ff = 1 / PixelSize * Zoom;
               xlim([ROI(1) - xExtra, ROI(2) + xExtra] * ff); 
               ylim([ROI(3) - yExtra, ROI(4) + yExtra] * ff); 
               if GaussianImageKludge
                  ylim([GaussianImageKludge - ROI(4) - yExtra, ...
                        GaussianImageKludge - ROI(3) + yExtra] * ff); 
               end
            else
               xlim([ROI(1) - xExtra, ROI(2) + xExtra]); 
               ylim([ROI(3) - yExtra, ROI(4) + yExtra]); 
               if GaussianImageKludge
                  ylim([GaussianImageKludge - ROI(4) - yExtra, ...
                        GaussianImageKludge - ROI(3) + yExtra]); 
               end
            end
         end

         title(sprintf('%s ROI %d', shrt, j));
         axis off
         hold off

         if opt.SR
            in_type = 'SR';
         elseif opt.BaGoL
            in_type = 'BaGoLSR';
         elseif opt.MAPN
            in_type = 'MAPN';
         elseif opt.Circle
            in_type = 'SR+MAPN';
         end

         if opt.Dot
            out_type = 'dot';
         elseif opt.Gaussian
            out_type = 'Gaussian';
         elseif opt.Circle
            out_type = 'circle';
         end

         if ~opt.NoSave
            SaveFile = fullfile(SaveDir, sprintf('%s_ROI%d_%s_%s', ...
                                              short, j, in_type, out_type));

            % Use imwrite rather than print for good resolution, but only for
            % the basic bitmap image.
            fprintf('\nSaving %s ...\n', SaveFile);
            if (opt.Gaussian || opt.Circle) && ~opt.Boundary && ~opt.Cluster
               imwrite(MapIm, [SaveFile, '.png']);
            else
               print(gcf, SaveFile, '-dpng', '-r600');
            end
            if opt.Dot
               saveas(gcf, SaveFile, 'fig');
            end
         else
            pause(3)
         end
         close
      end % for j (ROI #)
      close all
   end % for i (file #)

end
