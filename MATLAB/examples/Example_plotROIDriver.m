% A simple scipt to invoke the plotROIDriver (which invokes plotROI) in various
% ways.

close all

%plotROIDriver(PixelSize, options, start_datadir, SaveDir)

% One BaGoL MAPN image exhibiting multiple ROIs.
options = {'MAPN', 'Gaussian', 'Boundary', 'Cluster', 'OneImage', 'NoSave'};

% One BaGoL MAPN image per ROI.
%options = {'MAPN', 'Gaussian', 'Boundary', 'Cluster', 'ROIImages', 'NoSave'};
%options = {'MAPN', 'Gaussian', 'Cluster', 'ROIImages'};

% One BaGoL MAPN ROI per MAPN Results file.
%options = {'MAPN', 'MAPNResultsROI', 'Gaussian', 'Cluster'};
%options = {'MAPN', 'MAPNResultsROI', 'GaussSEConst', 'Cluster'};

% One SMITE SR image exhibiting multiple ROIs.
%options = {'SR', 'Gaussian', 'Boundary', 'Cluster', 'OneImage', 'NoSave'};
%options = {'SR', 'Dot', 'Boundary', 'Cluster', 'OneImage', 'NoSave'};

% One SMITE SR image per ROI.
%options = {'SR', 'Gaussian', 'Boundary', 'Cluster', 'ROIImages', 'NoSave'};
%options = {'SR', 'Gaussian', 'Cluster', 'ROIImages'};

start_DataDir = '/mnt/nas/cellpath/Genmab/Data/';

SaveDir = '/mnt/nas/lidkelab/Personal Folders/MJW/ROI/NEW/OUT';

% Cells to produce plots for using an absolute numbering over the cells that
% were selected,
IncludeCell = [];
%IncludeCell = [13, 25]; % DF3
%IncludeCell = [11, 20]; % DNP-BSA
%IncludeCell = [9, 116]; % DNP-PEG-BSA
%IncludeCell = [1, 8]; % DNP-BSA_10nM
%IncludeCell = [1, 7]; % DNP-BSA_2150pM
%IncludeCell = [2]; % Resting

PixelSize = 97.8;   % nm/pixel

% Plot dot, Gaussian or circle images of SMD/MAPN coordinates per ROI per Cell.
% In other words, this function can plot cluster boundaries overlaid on dot,
% Gaussian or circle plots of SR or BaGoL results.  The driver collects the
% paths and files needed, while plotROI produces the specified plots from these
% and a few other pieces of information for each ROI in each Cell, putting the
% results in SaveDir.  WARNING: Not ALL combinations of options work together!
%
% INPUTS:
%    PixelSize       plxel length (nm)
%    options         cell array of strings that are selected from the following
%                    options [DEFAULT = {'BaGoL', 'Gaussian', 'Cluster'}]:
%       SR              SR Results file
%       BaGoL           BaGoL Results file (BGL.SMD)
%       MAPN            BaGoL MAPN file
%       MAPNResultsROI  Individual BaGoL MAPN Results_ROI files
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

CI = smi_cluster.ClusterInterface;
CI.plotROIDriver(PixelSize, options, start_DataDir, SaveDir, IncludeCell);
