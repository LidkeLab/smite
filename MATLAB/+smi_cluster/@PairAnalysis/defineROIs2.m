function n_ROIs_ALL = ...
   defineROIs2(Files1, Files2, Pixel2nm, Color, ROI_sizes, ResultsDir, ...
               RegistrationNeeded, RegistrationViaDC)
% Select ROIs simultaneously for label 1 and label 2 over all images.
%
% INPUTS:
%    Files1       list of label 1 files fully specified
%    Files2       list of label 2 files fully specified
%                 NOTE: The number of label 1 and label 2 files should be the
%                 same as the ROIs will be specified for both simultaneously
%    Pixel2nm     conversion factor from pixels to nm
%    Color        label 1 and 2 colors on display
%    ROI_sizes    ROI length and width, e.g., [2000, 2000] (nm)
%    ResultsDir   directory in which to place the files generated
%    RegistrationNeeded   2-color registration via fiducials needed?
%    RegistrationViaDC    2-color registration via drift correction?
%
% OUTPUTS:
%    n_ROIs_ALL   total number of ROIs selected over all label 1/2 image pairs
%    *_ROIs.mat   data files specifying the ROIs selected for each image
%    *_ROIs.png   images identifying the ROIs selected per pair of label 1/2
%                 files
%    Note: data file defines n_ROIs (number of ROIs) and
%       RoI{1:n_ROIs}.ROI             ROI coordinates (xmin, xmax, ymin, ymax)
%       RoI{1:n_ROIs}.X/Y/X_SE/Y_SE   [L1 values, L2 values]
%       RoI{1:n_ROIs}.SMD             {L1 values, L2 values}

% Created by
%    Michael J. Wester (2022)

   RT = smi_helpers.ROITools();
   RT.EM = false;      % if true, data is in EM format
   RT.Color = Color;   % label 1 and 2 colors on display
   RT.ROI_sizes = ROI_sizes;
   RT.Pixel2nm = Pixel2nm;

   n_ROIs_ALL = 0;
   n_files = numel(Files1);
   for i = 1 : n_files
      Label1 = Files1{i};
      Label2 = Files2{i};

      B1 = load(Label1);
      B2 = load(Label2);
      if isfield(B1, 'BGL') || isfield(B1, 'MAPN')
         % BaGoL_Results files (coordinate units are nm).
         if isfield(B1, 'BGL')
            SMD1 = B1.BGL.MAPN;
            SMD2 = B2.BGL.MAPN;
         else
            SMD1 = B1.MAPN;
            SMD2 = B2.MAPN;
         end
         % Convert to pixels from nm as this is now an SMD structure.
         SMD1.X = SMD1.X / Pixel2nm;
         SMD1.Y = SMD1.Y / Pixel2nm;
         SMD2.X = SMD2.X / Pixel2nm;
         SMD2.Y = SMD2.Y / Pixel2nm;
         SMD1.X_SE = SMD1.X_SE / Pixel2nm;
         SMD1.Y_SE = SMD1.Y_SE / Pixel2nm;
         SMD2.X_SE = SMD2.X_SE / Pixel2nm;
         SMD2.Y_SE = SMD2.Y_SE / Pixel2nm;
      elseif isfield(B1, 'SMD')
         % SMD Results files (coordinate units are pixels).
         SMD1 = B1.SMD;
         SMD2 = B2.SMD;
      elseif isfield(B1, 'SMR')
         % SR ResultsStruct files.
         SMD1 = B1.SMR;
         SMD2 = B2.SMR;
      else
         error('Unknown result structure!');
      end

      if RegistrationNeeded
         [SMD2.X, SMD2.Y] = ...
            transformPointsInverse(R.LWMTransform, SMD2.X, SMD2.Y);
      end
      if RegistrationViaDC
         [delta12, Statistics] = ...
            smi_core.DriftCorrection.regViaDC(SMD1, SMD2);
         fprintf('regViaDC Cell %d deltas: [%g, %g]\n', i, delta12);
         SMD2.X = SMD2.X - delta12(1);
         SMD2.Y = SMD2.Y - delta12(2);
      end

      % The below assumes files are named like *_L1.ext and *_L2.ext where
      % L1/L2 are text distinguishing Label 1 from Label 2, and .ext is the
      % file extension (e.g., '.mat', 'txt', etc.).
      [filepath, filename] = fileparts(Files1{i});
      desc = regexprep(filename, '_ResultsStruct$', '');

      % Select ROIs for each pair of files (identified by Label1 and Label2).
      txt = regexprep(desc, '_', '\\_');
      %  Use the mouse to select ROIs (regions of interest):
      %    left click  chooses the center of a fixed size (x_size x y_size)
      %                region
      %    right click chooses an adjustable rectangular size region
      %    key press:
      %       backspace or delete   deletes the previous region
      %       anything else         terminates selection
      [n_ROIs, RoI, XYsize] = RT.getROI({SMD1, SMD2}, txt);
      n_ROIs_ALL = n_ROIs_ALL + n_ROIs;

      % Redefine Pixel2nm for later use.  This is done as original coordinates
      % are in pixels, but are now converted into nm for use in the RoI
      % structure, and uses such as clustering and related statistics.
      % This is only really needed for H-SET clustering.  (2024/12/11)
      Pixel2nmSAVE = Pixel2nm;
%     Pixel2nm = 1;
      saveas(gcf, fullfile(ResultsDir, sprintf('%s_ROIs.fig', desc)));
      print(fullfile(ResultsDir, sprintf('%s_ROIs.png', desc)), '-dpng');
      save(fullfile(ResultsDir, sprintf('%s_ROIs.mat', desc)), ...
           'Label1', 'Label2', 'Pixel2nm', 'XYsize', 'n_ROIs', 'RoI');
      close
      Pixel2nm = Pixel2nmSAVE;
      fid = fopen(fullfile(ResultsDir, sprintf('%s_ROIs.txt', desc)), 'w');
      fprintf(fid, 'ROI #, [xmin, xmax, ymin, ymax]\n');
      for i = 1 : n_ROIs
         fprintf(fid, '%d %7.3f %7.3f %7.3f %7.3f\n', ...
                      i, RoI{i}.ROI ./ Pixel2nm);
      end
      fclose(fid);
   end
   fprintf('Done ROIs.\n');

end
