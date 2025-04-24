function defineROIs(pathname, files, Pixel2nm, RT, oneROI, ROI_sizes)
% ----------- Define the ROIs
% Choose ROIs of a fixed size over a series of images.  These are typically
% used for cluster analysis.
%
% INPUTS:
%    pathname          path to where the files are located
%    files             cell array of _Results*.mat file names containing
%                      coordinate data
%    Pixel2nm          conversion factor from pixels to nm
%    RT                class reference to smi_helpers.ROITools()
%    oneROI            treat the whole image as a single ROI (logical)
%    ROI_sizes         ROI x, y sizes (in nm) used if oneROI is true
%
% OUTPUTS:
%    pathname/Analysis/*_ROIs.mat is saved containing the coordinates of the
%    ROIs selected for each cell image.  This selection process is performed
%    for each file name (cell image) provided, so the *_ROIs.mat files will be
%    identified by the names of the files containing the original (x, y)
%    coordinates.

% Created by
%    Michael J. Wester (2022)

   results_dir = fullfile(pathname, 'Analysis');
   % Create results_dir if it does not already exist.
   if ~isfolder(results_dir)
      mkdir(results_dir);
   end

   % Define the ROIs for each image.
   n_files = numel(files);
   for j = 1 : n_files
      short = regexprep(files{j}, '.mat$', '');
      ResultsFile = fullfile(pathname, files{j});
      fprintf('%s ...\n', short);
      if oneROI
         [XY, XY_STD, XYsize] = ...
            RT.import_XY(ResultsFile, Pixel2nm, '');
         n_ROIs = 1;
         RoI{1}.ROI = [0, ROI_sizes(1), 0, ROI_sizes(2)];
         RoI{1}.X = {XY(:, 1)};
         RoI{1}.Y = {XY(:, 2)};
         RoI{1}.X_STD = {XY_STD(:, 1)};
         RoI{1}.Y_STD = {XY_STD(:, 2)};
      else
         [n_ROIs, RoI, XYsize] = ...
            RT.getROI(ResultsFile, files{j});
         saveas(gcf, fullfile(results_dir, sprintf('%s_ROIs.fig', short)));
         saveas(gcf, fullfile(results_dir, sprintf('%s_ROIs.png', short)));
      end
      OriginLLvsUL = RT.OriginLLvsUL;
      save(fullfile(results_dir, [short, '_ROIs.mat']), ...
           'ResultsFile', 'Pixel2nm', 'XYsize', 'n_ROIs', 'RoI', ...
           'OriginLLvsUL');
      close
      fid = fopen(fullfile(results_dir, [short, '_ROIs.txt']), 'w');
      for i = 1 : n_ROIs
         fprintf(fid, '%d %7.3f %7.3f %7.3f %7.3f\n', ...
                      i, RoI{i}.ROI ./ Pixel2nm);
      end
      fclose(fid);
   end
   fprintf('Done ROIs.\n');

end
