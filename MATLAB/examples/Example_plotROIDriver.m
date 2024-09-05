% A simple scipt to invoke the plotROIDriver (which invokes plotROI) in various
% ways.

close all

%plotROIDriver(PixelSize, options, start_datadir, SaveDir)

% One BaGoL MAPN image exhibiting multiple ROIs.
options = {'MAPN', 'Gaussian', 'Boundary', 'Cluster', 'OneImage', 'NoSave'};

% One BaGoL MAPN image per ROI.
%options = {'MAPN', 'Gaussian', 'Boundary', 'Cluster', 'ROIImages', 'NoSave'};
%options = {'MAPN', 'Gaussian', 'Cluster', 'ROIImages'};

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

CI = smi_cluster.ClusterInterface;
CI.plotROIDriver(97.8, options, start_DataDir, SaveDir, IncludeCell);
