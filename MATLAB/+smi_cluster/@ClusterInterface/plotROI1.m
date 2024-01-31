function plotROI1(opt, pathnameC, filesC, pathnameB, filesB, PixelSize, ...
                  SaveDir)
% Plot dot, Gaussian or circle images of SMD/MAPN coordinates per ROI per cell.
% Just one image is produced per cell displaying all the ROIs.
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

      % Plot ROIs on top of one base image.
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

      % Moved out of the j loop.
      if opt.MAPN
         BGL.MAPN = SMDsave;
         if opt.Gaussian
            BGL.PixelSize = 1;
         end
         BGL.PImageSize = dataB.PImageSize;
         BGL.XStart = dataB.XStart;
         BGL.YStart = dataB.YStart;
         if GaussianImageKludge > 0
            BGL.MAPN.Y = GaussianImageKludge - BGL.MAPN.Y;
         end
         MapIm = CI.genMAPNIm1(BGL, 1);
         MapIm = smi.BaGoL.scaleIm(MapIm, 98);
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
%           if opt.MAPN
%              BGL.MAPN = SMD;
%           end
            % Gaussian image plot of SR localizations.
            BGL.PixelSize = 1;
            %if GaussianImageKludge > 0
            %   BGL.MAPN.Y = GaussianImageKludge - BGL.MAPN.Y;
            %end
            GIK = GaussianImageKludge;
            indx = ROI(1) - xExtra <= SMD.X & SMD.X <= ROI(2) + yExtra & ...
                   ROI(3) - yExtra <= GIK - SMD.Y &                      ...
                   GIK - SMD.Y <= ROI(4) + yExtra;
            SMD.X = SMD.X(indx);
            SMD.Y = SMD.Y(indx);
            SMD.X_SE = SMD.X_SE(indx);
            SMD.Y_SE = SMD.Y_SE(indx);
            if opt.BaGoL
               BGL.SMD = SMD;
               MapIm = CI.genMAPNIm1(BGL, 2);
%           elseif opt.MAPN
%              BGL.PImageSize = dataB.PImageSize;
%              BGL.XStart = dataB.XStart;
%              BGL.YStart = dataB.YStart;
%              if GaussianImageKludge > 0
%                 BGL.MAPN.Y = GaussianImageKludge - BGL.MAPN.Y;
%              end
%              MapIm = CI.genMAPNIm1(BGL, 1);
            end
            if ~opt.MAPN   % done above just once
               %MapIm = BGL.scaleIm(MapIm, 98);
               MapIm = smi.BaGoL.scaleIm(MapIm, 98);
            end
         end

         % Gaussian image plot of SR localizations.
         ScaleBarLength = 500;  % nm
         ScaleBarWidth  = 100;   % nm
         if opt.Gaussian
            if j == 1
               fprintf('ScaleBar = %d nm\n', ScaleBarLength);
               % Add in a scale bar (lower right).
               Xoffset = 500;
               Yoffset = 250;
               %Xstart = ROI(2) - Xoffset - ScaleBarLength;
               %Ystart = ...
               %   GaussianImageKludge - (ROI(3) + Yoffset + ScaleBarWidth);
               Xstart = dataB.PImageSize - Xoffset - ScaleBarLength;
               Ystart = ...
                  GaussianImageKludge - (0 + Yoffset + ScaleBarWidth);
               X = round(Xstart) : round(Xstart + ScaleBarLength);
               Y = round(Ystart) : round(Ystart + ScaleBarWidth);
               MapIm(Y, X) = 255;

               fprintf('ROI (%d):', dataC.n_ROIs(i));
               % MATLAB is acting weird here.  I can display interactively
               % multiple figures with the code below, but it screws up when
               % saving the images (hence the opt.NoSave).  However, saving the
               % images correctly does not produce multiple interactive figures
               % (just the last one).  So one has choose interactive figures vs
               % saved images (which don't look that great anyway because I
               % cannot save an image overlayed with points and lines at a high
               % resolution, however, it displays just fine.  In short, images
               % and plot figures don't play well together!
               if i > 1 && opt.NoSave
                  figure
               end
               imshow(MapIm, hot);
               hold on
            end
            fprintf(' %d', j);
            if j == dataC.n_ROIs(i)
               fprintf('\n');
            end
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
            plot(x, y, 'g-', 'LineWidth', 2);
         end

         if  yExtra > 0
            % Add a label designating the ROI number.
            x_label = (ROI(1) + ROI(2))/2;
	    %y_label = ROI(4) + yExtra/2;
            y_label = (ROI(3) + ROI(4))/2;
            if GaussianImageKludge > 0
               y_label = GaussianImageKludge - y_label;
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
               C = results.C;
               if numel(C{l}) <= 2
                  %plot(x(C{l}), y(C{l}), 'm.-', 'LineWidth', 2, ...
                  plot(x(C{l}), y(C{l}), 'm-', 'LineWidth', 1, ...
                       'MarkerSize', 12);
               else
                  % Determine cluster boundary indices from the cluster
                  % contents.
                  k = boundary(double(x(C{l})), double(y(C{l})), ShrinkFactor);
                  k = C{l}(k);
                  k = [k, k(1)];
                  %plot(x(k), y(k), 'm.-', 'LineWidth', 2, 'MarkerSize', 12);
                  plot(x(k), y(k), 'm-', 'LineWidth', 1, 'MarkerSize', 12);
               end
            end % for l (cluster #)
         end % if opt.Cluster

% Axis limits no longer needed for the full image.
%        % Circle plots already are the full ROI.
%        if ~opt.Circle
%           xlim([ROI(1) - xExtra, ROI(2) + xExtra]); 
%           ylim([ROI(3) - yExtra, ROI(4) + yExtra]); 
%           if GaussianImageKludge
%              ylim([GaussianImageKludge - ROI(4) - yExtra, ...
%                    GaussianImageKludge - ROI(3) + yExtra]); 
%           end
%        end

%        title(sprintf('%s ROI %d', shrt, j));
%        axis off
%        hold off

%        if opt.SR
%           in_type = 'SR';
%        elseif opt.BaGoL
%           in_type = 'BaGoLSR';
%        elseif opt.MAPN
%           in_type = 'MAPN';
%        elseif opt.Circle
%           in_type = 'SR+MAPN';
%        end

%        if opt.Dot
%           out_type = 'dot';
%        elseif opt.Gaussian
%           out_type = 'Gaussian';
%        elseif opt.Circle
%           out_type = 'circle';
%        end

%        SaveFile = fullfile(SaveDir, sprintf('%s_ROI%d_%s_%s', ...
%                                             short, j, in_type, out_type));

%        % Use imwrite rather than print for good resolution, but only for the
%        % basic bitmap image.
%        if (opt.Gaussian || opt.Circle) && ~opt.Boundary && ~opt.Cluster
%           imwrite(MapIm, [SaveFile, '.png']);
%        else
%           print(gcf, SaveFile, '-dpng', '-r600');
%        end
%        if opt.Dot
%           saveas(gcf, SaveFile, 'fig');
%        end
%        close
      end % for j (ROI #)

      title(shrt);
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
         SaveFile = fullfile(SaveDir, sprintf('%s_%s_%s', ...
                                              short, in_type, out_type));
         fprintf('Saving %s ...\n', SaveFile);
         print(gcf, SaveFile, '-dpng', '-noui', '-r1200');
         %cdata = print(gcf, '-RGBImage', '-noui', '-r1200');
         %imwrite(cdata, [SaveFile, '.png']);
         %print(gcf, SaveFile, '-dpdf', '-noui', '-r600');
      else
         pause(3)
      end

%     close all
   end % for i (file #)

end
