function overlayBaGoLROIs(pathnameB, filesB, MAPNfile, ROI_sizes, SRImageZoom)
% ---------- Possibly, produce Gaussian overlay images of individually
%            processed BaGoL ROIs
% This function will overlay label 1 onto label 2 Gaussian images produced from
% the input ROI file names into appropriately named files.  This assumes the
% various files involved are following a naming convention.  For example, the
% BaGoL results files oer ROI are named:
%    BaGoL_Results_Cell_01_Label_01_Results_ROI_01_ResultsStruct.mat
%    BaGoL_Results_Cell_01_Label_01_Results_ROI_02_ResultsStruct.mat
%    ...
%    BaGoL_Results_Cell_01_Label_02_Results_ROI_01_ResultsStruct.mat
%    ...
%    BaGoL_Results_Cell_02_Label_01_Results_ROI_01_ResultsStruct.mat
%    ...
% and the BaGoL MAPN coordinates per ROI are named:
%    MAPN_Cell_01_Label_01_Results_ROI_01.mat
%    MAPN_Cell_01_Label_01_Results_ROI_02.mat
%    ...
%    MAPN_Cell_01_Label_02_Results_ROI_01.mat
%    ...
%    MAPN_Cell_02_Label_01_Results_ROI_01.mat
%    ...
% The collection of ROIs (overlaying both labels) produced from the BaGoL data:
%    Cell_01_ROI_01_OverlayLabels_BaGoL.png
%    Cell_01_ROI_02_OverlayLabels_BaGoL.png
%    ...
%    Cell_02_ROI_01_OverlayLabels_BaGoL.png
%    ...
%
% INPUTS:
%    pathnameB   path to where the BaGoL files below are located
%    filesB      BaGoL files to collect ROI coordinate info from
%    MAPNfile    if false, assume BaGoL_Results_*_Results*.mat files,
%                otherwise if true, assume MAPN_*.mat files
%    ROI_sizes   ROI dimensions in x and y (nm)
%    SRImageZoom image resolution magnification factor [default = 4]
%
% OUTPUTS:
%    Saves pathnameB/Analysis/*_OverlayLabels_BaGoL.png
%    Saves pathnameB/Analysis/*_OverlayLabels_BaGoL_circle.png

% Created by
%    Michael J. Wester (2023)

   if ~exist('SRImageZoom', 'var')
      SRImageZoom = 4;
   end

   results_dir = pathnameB;

   n_labels = 2;
   n_Label2 = 0;
   n_files = numel(filesB);
   % Each file contains a singly labeled ROI.
   for i = 1 : n_files
      fileB = filesB{i};
      % Extract the cell, label and ROI numbers from the file name.
      c = str2num(regexprep(fileB, '^.*Cell_([0-9][0-9]).*$',  '$1'));
      l = str2num(regexprep(fileB, '^.*Label_([0-9][0-9]).*$', '$1'));
      r = str2num(regexprep(fileB, '^.*ROI_([0-9][0-9]).*$',   '$1'));
      fprintf('Cell %02d ROI %02d Label %02d\n', c, r, l);

      % If this is a label 2 file name, assume it has been or will be processed
      % with the label 1 file, so skip this iteration of the loop.
      if l == 2
         n_Label2 = n_Label2 + 1;
         continue;
      end

      ResultsFile = sprintf('Cell_%02d_ROI_%02d_OverlayLabels_BaGoL', c, r);
      SaveFile = fullfile(pathnameB, ResultsFile);
      ResultsFile = sprintf('Cell_%02d_ROI_%02d_OverlayLabels_BaGoL_circle',...
                            c, r);
      SaveFileC = fullfile(pathnameB, ResultsFile);

      % Here, fileB is label 1.
      fileB1 = fileB;
      fileB2 = regexprep(fileB1, 'Label_01', 'Label_02');

      dataB1 = load(fullfile(pathnameB, fileB1));
      dataB2 = load(fullfile(pathnameB, fileB2));

      ScalebarLength = 500; % nm (smi.BaGoL.scalebar)
      if MAPNfile
         dataB1_MAPN = dataB1.MAPN;
         dataB2_MAPN = dataB2.MAPN;
         XStart = dataB1.XStart;
         YStart = dataB1.YStart;
         PImageSize = dataB1.PImageSize;
         PixelSize = dataB1.PixelSize;
         PixelSize1 = PixelSize;

         BGL1 = smi.BaGoL;
         BGL1.XStart = XStart;
         BGL1.YStart = YStart;
         BGL1.PImageSize = PImageSize;
         BGL1.PixelSize = PixelSize / SRImageZoom;
         BGL1.MAPN = dataB1_MAPN;

         BGL2 = smi.BaGoL;
         BGL2.XStart = XStart;
         BGL2.YStart = YStart;
         BGL2.PImageSize = PImageSize;
         BGL2.PixelSize = PixelSize / SRImageZoom;
         BGL2.MAPN = dataB2_MAPN;
      else
         dataB1_MAPN = dataB1.BGL.MAPN;
         dataB2_MAPN = dataB2.BGL.MAPN;
         XStart = dataB1.BGL.XStart;
         YStart = dataB1.BGL.YStart;
         PixelSize = dataB1.BGL.PixelSize;

         BGL1 = dataB1.BGL;
         PixelSize1 = BGL1.PixelSize;
         BGL1.PixelSize = BGL1.PixelSize / SRImageZoom;

         BGL2 = dataB2.BGL;
         PixelSize2 = BGL2.PixelSize;
         BGL2.PixelSize = BGL2.PixelSize / SRImageZoom;
      end

      % Produce Gaussian Images for the two labels, then overlay them.
      BaGoLScalebarLength = ScalebarLength / BGL1.PixelSize;
      GaussIm1 = BGL1.genMAPNIm(1);
      GaussIm1 = BGL1.scaleIm(GaussIm1, 98);
      GaussIm1 = BGL1.scalebar(GaussIm1, 1, BaGoLScalebarLength);

      BaGoLScalebarLength = ScalebarLength / BGL2.PixelSize;
      GaussIm2 = BGL2.genMAPNIm(1);
      GaussIm2 = BGL2.scaleIm(GaussIm2, 98);
      GaussIm2 = BGL2.scalebar(GaussIm2, 1, BaGoLScalebarLength);

      ImageStack = zeros([size(GaussIm1, 1), size(GaussIm1, 2), 2]);
      ImageStack(:, :, 1) = GaussIm1;
      ImageStack(:, :, 2) = GaussIm2;
      [OverlayImage, ColorOrderTag] = ...
         smi_vis.GenerateImages.overlayNImages(ImageStack);
      imwrite(OverlayImage, hot(256), [SaveFile, '.png']);

      % Circle image of label 1 versus label 2.
      XSize = ROI_sizes(1);
      YSize = ROI_sizes(2);
      PixelSize = PixelSize1 / SRImageZoom;
      RadiusScale = 2;
      CircleIm = BGL1.genSRMAPNOverlay(BGL1.MAPN, BGL2.MAPN, ...
                    XSize, YSize, PixelSize, pathnameB,      ...
                    BGL1.XStart, BGL1.YStart, RadiusScale, ScalebarLength);
      imwrite(CircleIm, [SaveFileC, '.png']);
   end

   % Consistency check.
   if n_labels * n_Label2 ~= n_files
      error('%d (# labels) * %d (# label 2 ROIs) != %d (# files)!', ...
            n_labels, n_Label2, n_files);
   end

   fprintf('Done overlaying BaGoL ROIs: %s.\n', results_dir);

end
