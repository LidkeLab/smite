function plotROI(opt, pathnameC, filesC, pathnameB, filesB, PixelSize, SaveDir)
% Plot dot, Gaussian or circle images of SMD/MAPN coordinates per ROI per cell.
% Special case when individual BaGoL MAPN Results_ROI files were created
% (MAPNResultsROI = true).  These are special files, so this routine (derived
% from plotROI.m) has been simplified; the x'ed out options below are not valid
% here and will be ignored.  Circle plots only have MAPN coordinates available
% to them for now, which are colored green.
%
% Each individual (cell #, ROI #) pair in filesB is processed one at a time,
% making one of the plot types (Gaussian, GaussSEConst, Circle), optionally
% overlayed by cluster boundaries.  Images are saved unless NoSave is
% specified.
%
% filesB and filesC.files are assumed to have the string 'Cell_nn' embedded in
% their names, where nn is a 2-digit cell number like 01, 12, etc.  filesB is
% also assumed to have a 2-digit ROI number: 'ROI_nn' embedded in their names.
%
% INPUTS:
%    opt         data characteristics and types of plots to produce if true
% x     SR             SR Results file
% x     BaGoL          BaGoL Results file (BGL.SMD)
%       MAPN           BaGoL MAPN file
%       MAPNResultsROI Individual BaGoL MAPN Results_ROI files
% x     Dot            Dot plot
%       Gaussian       Gaussian plot
%       GaussSEConst   Gaussian plot with constant X/Y_SE
%       Circle         Circle plot (BaGoL Results file: BGL.SMD + BGL.MAPN)
% x     Boundary       Include ROI boundaries
%       Cluster        Include ROI clusters
%       NoSave         Do not save outputs
%
%       IncludeCell    Cells to produce plots for [DEFAULT: 1 : numel(filesB)]
%                      using the numbering in the list filesB
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
%    Michael J. Wester (2025)

   % If zero, ROIs in simpleROIcluster were chosen using plotted points, but
   % if > 0, ROIs were chosen using GaussianIm = true so following the DIPimage
   % convention---the value should then be YSize * PixelSize.
   % Changes (9/20/21) in simpleROIcluster seem to have made this fix obsolete,
   % but will retain the code for new just in case they are needed once again.
   % Certainly useful when plotting boundaries/clusters on top of images.
   %GaussianImageKludge = 256*PixelSize;
   % Boundary shrink factor (0 = convex hull, 1 = as concave as possible,
   % 0.5 = MATLAB default)
   ShrinkFactor = 0.5;
   ScaleBarLength = 500; % nm
   ScaleBarWidth  = 100; % nm
   % Special BaGoL plots for Diane: X/Y_SE are constant and uniform.
   SEConstant = 5; % nm

   CI = smi_cluster.ClusterInterface();

   dataC = load(fullfile(pathnameC, filesC{1}));

   ScaleBarMsg = true;

   % Extract the Cell numbers from the dataC filenames for matching purposes
   % with dataB Cell numbers (comes in when there are missing Cells.
   n_dataC = numel(dataC.files);
   c_CellNum = zeros(1, n_dataC);
   c_CellNum2Index = zeros(1, n_dataC);
   for i = 1 : n_dataC
      c_CellNum(i) = ...
         str2num(regexprep(dataC.files{i}, '^.*Cell_([0-9]+).*$', '$1'));
      c_CellNum2Index(c_CellNum(i)) = i;
   end

   for i = 1 : numel(filesB)
      if ~ismember(i, opt.IncludeCell)
         % Remember: j = sum(dataC.n_ROIs(1 : i_cell - 1)) + i_ROI;
         % Kind of kludgy, but consolidated ROIs use ALL the cells in filesB to
         % determine their numbering.
         continue;
      end

      fileB = filesB{i};
      i_cell = str2num(regexprep(fileB, '^.*Cell_([0-9]+).*$', '$1'));
      i_ROI  = str2num(regexprep(fileB, '^.*ROI_([0-9]+).*$',  '$1'));
      dataB = load(fullfile(pathnameB, fileB));
%     SMD  = dataB.SMD;
      SMD  = dataB.MAPN;
      MAPN = dataB.MAPN;
      % Remove extraneous material from the filenames.
      short = regexprep(fileB, '.mat$', '');
      short = regexprep(short, '_ResultsStruct$', '');
      short = regexprep(short, '_Results$', '');
      short = regexprep(short, '^BaGoL_Results_', '');
      shrt  = regexprep(short, '_', '\\_');

      % Display a quick summary of the upcoming analysis.
      if opt.MAPN
         fprintf('%s: MAPN = %d\n', short, numel(MAPN.X));
      end

      % Plot ROIs individually.
      if opt.MAPN
         SMDsave  = MAPN;
         MAPNsave = MAPN;
      elseif opt.Circle
         SMDsave  = BGL.SMD;
         MAPNsave = BGL.MAPN;
      end
      if opt.GaussSEConst
         SMDsave.X_SE = single(SEConstant) .* ones(size(SMDsave.X));
         SMDsave.Y_SE = SMDsave.X_SE;
      end

      % Added indexing arrays to deal wih missing Cells.
      %j = sum(dataC.n_ROIs(1 : i_cell - 1)) + i_ROI;
      %ROI = dataC.RoI{i_cell}{i_ROI}.ROI;
      j = 0;
      for k = 1 : i_cell - 1
         if 0 < c_CellNum(k) & c_CellNum(k) < i_cell
            j = j + dataC.n_ROIs(k);
         end
      end
      j = j + i_ROI;
      ROI = dataC.RoI{c_CellNum2Index(i_cell)}{i_ROI}.ROI;
      SMD = SMDsave;

      if opt.Gaussian
         if opt.MAPN
            BGL.MAPN = SMD;
         end
         % Gaussian image plot of SR localizations.
         BGL.PixelSize = 1;
         if opt.MAPN
            BGL.PImageSize = dataB.PImageSize;
            BGL.XStart = dataB.XStart;
            BGL.YStart = dataB.YStart;
            MapIm = CI.genMAPNIm1(BGL, 1);
         end
         MapIm = smi.BaGoL.scaleIm(MapIm, 98);
      end

      fprintf('Cell %d ROI %d:', i_cell, i_ROI);
      fprintf(' (abs ROI %d)', j);

      % Gaussian image plot of SR localizations.
      ScaleBarLength = 100;  % nm
      ScaleBarWidth  = 25;   % nm
      if opt.Gaussian
         if ScaleBarMsg
            fprintf('\nScaleBar = %d nm', ScaleBarLength);
            ScaleBarMsg = false;
         end

         % Add in a scale bar (lower right).
         Xoffset = 100;
         Yoffset = 100;
         Xstart = ROI(2) - ROI(1) - (Xoffset + ScaleBarLength);
         Ystart = ROI(4) - ROI(3) - (Yoffset + ScaleBarWidth);
         X = round(Xstart) : round(Xstart + ScaleBarLength);
         Y = round(Ystart) : round(Ystart + ScaleBarWidth);
         MapIm(Y, X) = 255;
         imshow(MapIm, hot);

         hold on
      elseif opt.Circle
         BGL.PixelSize = 1;
         BGL.XStart = dataB.XStart;
         BGL.YStart = dataB.YStart;
         SMD = SMDsave;
         BGL.SMD = SMD;

         MAPN = MAPNsave;
         BGL.MAPN = MAPN;

         ROISize = ROI(2) - ROI(1);
         MapIm = CI.genSRMAPNOverlay1(BGL.SMD, BGL.MAPN, ROISize,      ...
                                      ROISize, 1, SaveDir, BGL.XStart, ...
                                      BGL.YStart, 1, ScaleBarLength);

         % Gaussian image plot of SR (green) + MAPN (magenta) localizations.
         imshow(MapIm, hot);
         hold on
      end

      if opt.Cluster
         % Plot cluster boundaries.
         results = dataC.results{j};
         for l = 1 : results.nC
            x = results.XY(:, 1) - BGL.XStart;
            y = results.XY(:, 2) - BGL.YStart;
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

      %title(sprintf('%s ROI %d', shrt, j));
      title(sprintf('%s', shrt));
      axis off
      hold off

      if opt.MAPN
         in_type = 'MAPN';
      elseif opt.Circle
         in_type = 'SR+MAPN';
      end

      if opt.Gaussian
         out_type = 'Gaussian';
         if opt.GaussSEConst
            out_type = 'GaussSEConst';
         end
      elseif opt.Circle
         out_type = 'circle';
      end

      if ~opt.NoSave
         %SaveFile = fullfile(SaveDir, sprintf('%s_ROI%d_%s_%s', ...
         %                                 short, j, in_type, out_type));
         SaveFile = fullfile(SaveDir, sprintf('%s_%s_%s', ...
                                           short, in_type, out_type));

         % Use imwrite rather than print for good resolution, but only for
         % the basic bitmap image.
         fprintf('\nSaving %s ...\n', SaveFile);
         if (opt.Gaussian || opt.Circle) && ~opt.Cluster
            imwrite(MapIm, hot, [SaveFile, '.png']);
         else
            print(gcf, SaveFile, '-dpng', '-r600');
         end
         if opt.Dot
            saveas(gcf, SaveFile, 'fig');
         end
      else
         fprintf('\n');
         pause(3)
      end
      close
   
   end % for i (file #)

end
