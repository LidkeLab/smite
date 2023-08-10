% Perform cluster analysis for comparison of experimental conditions.  This is
% a batch version of simpleROIcluster_smite for performing parameter studies.
% See pathname/files definition below in Statistics for a single condition.
% 
% NOTE: smite and PlotSpread should be on your MATLAB path.
%
% Work flow: (note results are put in the subdirectory 'Analysis')
%   Set important parameters:
%      ROI_sizes, Pixel2nm, start_datadir, oneROI
%   Statistics for a single condition:
%            _ROIs.mat -> ALL_*
%      (compute statistics for a single condition consisting of a number of
%       images all located in the same directory)
%
% Avoid spaces in filenames and conditions!!!

%% ---------- Set important parameters

doHopkins = true;  % Hopkins' test can be time consuming for dense ROIs 
% Set to false for very dense ROIs to avoid crashes due to lack of memory.
doSigmaActual = false;

ROI_sizes = [2000, 2000];   % [delta_x, delta_y] (nm)
A_ROI = prod(ROI_sizes);    % ROI area (nm^2)
%Pixel2nm = 16000/150;       % conversion factor from pixels to nm
%Pixel2nm = 108.018;         % pixels to nm [TIRF]
Pixel2nm = 97.8;            % pixels to nm [sequential]

% Often, for BaGoL analyses, it is simpler to use a single, large, encompassing
% ROI rather than a series of small ROIs.
oneROI = false;
%oneROI = true;
if oneROI
   ROI_sizes = [256, 256] * Pixel2nm;   % (nm)
   A_ROI = prod(ROI_sizes);
end

% Select the files starting from start_datadir.
%start_datadir = 'NEEDS_TO_BE_SET!';
start_datadir = '.';

fprintf('Done set parameters.\n');

%  - - - - - Choose _ROIs.mat files

% Interactive input.
%[pathname, files] = smi_helpers.selectFiles(start_datadir, ...
%                       '_ROIs.mat files', '*_ROIs.mat');
%answer = inputdlg('Name of combined results file:');
%condition = answer{1};

% Look for condition/Analysis/*_ROIs.mat or uncomment "files = ..."
% lines below and list files manually.
% Use specified condition.
%condition = '2F8';
%condition = 'B12';
%condition = 'E345R';
condition = 'RGY';
pathname = fullfile(start_datadir, condition, 'Analysis');
FILES = dir(fullfile(pathname, '*_ROIs.mat'));
files = { FILES.name };

% Manual specification.
%files = {
%'Cell_01_Label_01_Results_ROIs.mat'
%};

% base_name is a short identifier describing the analysis.
base_name = condition;

%  - - - - - Set up batch analysis ranges

% Minimum distance between clusters or maximum distance between points in a
% cluster.
E_range = [30, 40, 50];   % (nm)
% Minimum number of points in a cluster.
minPts_range = [3, 6, 10];
% Clustering algorithm.
%algorithm_range = {'DBSCAN', 'Hierarchical', 'Voronoi'};
algorithm_range = {'DBSCAN'};

CI = smi_cluster.ClusterInterface();
CI.singleCondition(pathname, files, algorithm_range, E_range, minPts_range,...
                   Pixel2nm, base_name, A_ROI, doHopkins, doSigmaActual);
