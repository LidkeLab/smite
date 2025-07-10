function plotROIDriver(PixelSize, options, start_datadir, SaveDir, IncludeCell)
% Plot dot, Gaussian or circle images of SMD/MAPN coordinates per ROI per Cell.
% In other words, this function can plot cluster boundaries overlaid on dot,
% Gaussian or circle plots of SR or BaGoL results.  The driver collects the
% paths and files needed, while plotROI produces the specified plots from these
% and a few other pieces of information for each ROI in each Cell, putting the
% results in SaveDir.
%
% This function is now incorporated in ClusterInterface.
%
% INPUTS:
%    PixelSize       plxel length (nm)
%    options         cell array of strings that are selected from the following
%                    options [DEFAULT = {'BaGoL', 'Gaussian', 'Cluster'}]:
%       SR              SR Results file
%       BaGoL           BaGoL Results file (BGL.SMD)
%       MAPN            BaGoL MAPN file
%       Dot             Dot plot
%       Gaussian        Gaussian plot
%       GaussSEConst    Gaussian plot with constant X/Y_SE
%       Circle          Circle plot (BaGoL Results file: BGL.SMD + BGL.MAPN)
%       Boundary        Include ROI boundaries
%       Cluster         Include ROI clusters
%       ROIImages       produce one image per ROI
%       OneImage        produce one image displaying all ROIs
%       NoSave          do not save outputs
%    start_datadir   starting data directory from which SaveDir will be made a
%                    subdirectory of by default if SaveDir not provided
%                    [DEFAULT = '.']
%    SaveDir         directory in which to save images
%                    [DEFAULT = fullfile(start_datadir, 'ROIClusterAnalysis')]
%    IncludeCell     Cells to produce plots for [DEFAULT: 1 : numel(filesB)]
%                    using the numbering in the list filesB

   if ~exist('options', 'var')
      options = {'MAPN', 'Gaussian', 'Cluster', 'ROIImages'};
   end

   opt.SR = false;
   opt.BaGoL = false;
   opt.MAPN = false;
   opt.Dot = false;
   opt.Gaussian = false;
   opt.GaussSEConst = false;
   opt.Circle = false;
   opt.Boundary = false;
   opt.Cluster = false;
   opt.ROIImages = false;
   opt.OneImage = false;
   opt.NoSave = false;

   if contains('SR', options)
      opt.SR = true;
   end
   if contains('BaGoL', options)
      opt.BaGoL = true;
   end
   if contains('MAPN', options)
      opt.MAPN = true;
   end
   if contains('Dot', options)
      opt.Dot = true;
   end
   if contains('Gaussian', options)
      opt.Gaussian = true;
   end
   if contains('GaussSEConst', options)
      opt.GaussSEConst = true;
      opt.Gaussian = true;
   end
   if contains('Circle', options)
      opt.Circle = true;
   end
   if contains('Boundary', options)
      opt.Boundary = true;
   end
   if contains('Cluster', options)
      opt.Cluster = true;
   end
   if contains('ROIImages', options)
      opt.ROIImages = true;
   end
   if contains('OneImage', options)
      opt.OneImage = true;
   end
   if contains('NoSave', options)
      opt.NoSave = true;
   end

   CI = smi_cluster.ClusterInterface();

   if ~exist('start_datadir', 'var')
      start_datadir = '.';
   end
   if ~exist('SaveDir', 'var')
      SaveDir = fullfile(start_datadir, 'ROIClusterAnalysis');
   end
   if ~isfolder(SaveDir)
      mkdir(SaveDir);
   end

   % Clusters collected together in a single *_results.mat file.
   [pathnameC, filesC] = smi_helpers.selectFiles(start_datadir, ...
      '*_results.mat single file', '*_results.mat');
   fprintf('C: %s\n', fullfile(pathnameC, filesC{1}));
   if opt.SR
      [pathnameB, filesB] = smi_helpers.selectFiles(start_datadir, ...
         'SR _Results.mat files', '*_Results.mat');
   elseif opt.BaGoL 
      % BaGoL Result/ResultStruct files.
      [pathnameB, filesB] = smi_helpers.selectFiles(start_datadir, ...
         'MF BaGoL _Results*.mat files', 'BaGoL_Results_*_Results*.mat');
   elseif opt.MAPN
      % BaGoL MAPN files.
      [pathnameB, filesB] = smi_helpers.selectFiles(start_datadir, ...
         'MAPN*.mat files', 'MAPN_*.mat');
   end
   fprintf('B: %s\n', pathnameB);

   % Include all the ROIs in only some cells.
   if exist('IncludeCell', 'var') && ~isempty(IncludeCell)
      opt.IncludeCell = IncludeCell;
   else
      opt.IncludeCell = 1 : numel(filesB);
   end

   if opt.OneImage
      CI.plotROI1(opt, pathnameC, filesC, pathnameB, filesB, PixelSize, ...
                  SaveDir);
   end
   if opt.ROIImages
      CI.plotROI(opt, pathnameC, filesC, pathnameB, filesB, PixelSize, ...
                 SaveDir);
   end

   fprintf('Done.\n');

end
