% Work flow: (note results are put in the subdirectory 'Analysis')
% NOTE: smite should be on your MATLAB path.
%
%   Set important parameters:
%      ROI_sizes, Pixel2nm, start_datadir
%   Define the ROIs:
%            _ResultsStruct.mat (or _Results.mat) -> _ROIs.mat
%      (each image [ResultsStruct or Results file] produces a separate ROIs
%       file containing all the ROIs for that image; the idea is to do all the
%       ROI selection early on so it need not be repeated---each image can have
%       multiple ROIs; note that this deals with files in a single directory)
%   Analyze the ROIs one-by-one:
%            _ROIs.mat -> _results.mat, etc.
%      (compute results for each ROI defined by the selected files located in
%       the same directory)
%   Analyze a group of files of ROIs all together:
%            multiple _ROIs.mat
%                 -> ALL_ and analysis_ALL_*.mat
%      (compute combined results for the ROIs defined by the selected files
%       located in the same directory)
%
% Avoid spaces in filenames and conditions!!!
%
% See doAnalysis    for analyses to perform:
%     doBiStats     for mutual pairwise distance and bivariate Ripley's plots
%     doClustering  for cluster parameters: E and minPts
%     doClusterSep2 for cluster separation plots
%     doLocSep2     for nearest neighbor localization separation plots
%     doOverlap1    for overlaps between label 1 clusters and label 2 locs
%     doOverlap2    for overlaps between label 2 clusters and label 1 locs
%     doPairCorr    for pair correlation settings
%     doPlot2       for 2-color plots
%
% Output figure naming conventions:
%    *_ROI*_L1/2_algorithm            clustering by algorithm
%    *_ROI*_L1/2_pairwiseCDF/PDF      NN distances compared to a random dist.
%    *_ROI*_L1,L2_pairwisePDF2/CDF2   NN distances 2-label PDF/CDF
%    *_ROI*_bivripley                 bivariate Ripley's
%    *_ROI*_crosscorr                 ROIwise pairwise cross-correlations
%    *_RC_crosscorrR                  ROI combined pairwise cross-correlation
%    *_cL1_lL2                        L1 clusters, L2 localizations overlap
%    *_cL2_lL1                        L2 clusters, L1 localizations overlap
%    *_ROI*_L1+L2                     2-label plots of localizations
%    *_cs_RC_PDF/CDF                  ROI combined 2-label cluster separations
%    *_ls_ROI*                        histogram per ROI of loc. L1+L2 seps.
%    *_ls_RC_PDF/CDF                  ROI combined localization L1+L2 seps.

%% --- Set important parameters
% Collect the set of files expressing each label together.

PA = smi_cluster.PairAnalysis();

%start_datadir = '/mnt/nas/cellpath/Genmab/Data/20221119_CHO_HAEGFR_antiHA_antiC1Q/PairwiseCorr';
start_datadir = '.';
Analysis = 'Analysis';

particles = {'L1', 'L2'};     % label short descriptors
Color = ['g', 'm'];           % label 1 and 2 colors on display
%Pixel2nm = 108.018;           % pixels to nm [TIRF]
Pixel2nm = 97.8;              % pixels to nm [sequential]
HistBinSize = 5;              % # pixels per bin to collect correlation stats
%RmaxAxisLimit = -1;           % sets r axis limit for pair correlation plots
RmaxAxisLimit = 400;          % sets r axis limit for pair correlation plots
                              % if > 0 (nm)
ROI_sizes = [7000, 7000];     % ROI sizes (nm)
RegistrationNeeded = false;   % 2-color registration via fiducials needed?
RegistrationViaDC = false;    % 2-color registration via drift correction?
E = [20, 20];                 % epsilon: max distance between 2 adjacent points
                              % in a cluster per label (nm)
minPts = [3, 3];              % minimum number of points in a cluster per label
%property.Fig_ext = '';       % if blank, produce .fig files rather than .png

% ROIs to exclude for combined processing---for example for files 1, 2, 3
%    ExcludeROI = {{2, 1}, {3, [2, 4]}}
% excludes ROI 1 from file2, ROIs 2+4 from file3
ExcludeROI = {};

% Preprocess exclusions.
excludedFiles = [];
if ~isempty(ExcludeROI)
   excludedFiles = arrayfun(@(k) ExcludeROI{k}{1}, 1:numel(ExcludeROI));
end

% Set up the results directory.
results_top = uigetdir(start_datadir, 'Directory to hold Analysis directory');
results_dir = fullfile(results_top, Analysis);
% Create results_dir if it does not already exist.
if ~isdir(results_dir)
   mkdir(results_dir);
end
%data_dir = ...
%   uigetdir(start_datadir, 'Data directory (containing *_ROis.mat files)');
%data_dir = results_top;
data_dir = results_dir;
property.Results = results_dir;

% If true, look for MAPN_*.mat, otherwise *_Results*.mat for BaGoL coordinates.
MAPNfile = false;

% For overlap plots, do NOT plot non-overlapping localizations if false.
PlotNonOverlap = false;

% options:
%    'combined'    only produce combined results over all ROIs
%    'plotting'    produce most individual plots
%    'BiStats'     pairwise mutual distances and bivariate Ripley's statistics
%                  for each ROI
%    'Clustering'  clusters for each label per ROI
%    'Clustering2' C2C nearest neighbor distances between label 1/label 2
%                  clusters per ROI and combined over all ROIs
%    'LocSep2'     nearest neighbor localization distances between labels 1 & 2
%    'Overlap1'    overlaps between label 1 clusters and label 2 localizations
%    'Overlap2'    overlaps between label 2 clusters and label 1 localizations
%    'PairCorr'    pair correlation per ROI and combined over all ROIs
%    'Plot2'       2D plot per ROI
options = ["BiStats", "Clustering", "Clustering2", "LocSep2", ...
           "Overlap1", "Overlap2", "PairCorr", "Plot2"];
%options = ["BiStats", "Clustering", "LocSep2", "Overlap2", "PairCorr"];
%options = ["PairCorr"];

fprintf('Done set parameters.\n');

%% ---------- Define the ROIs

% ROIS in a cell are defined by n_ROIs (number of ROIs) and
%    RoI{1:n_ROIs}.ROI             ROI coordinates (xmin, xmax, ymin, ymax)
%    RoI{1:n_ROIs}.X/Y/X_SE/Y_SE   [L1 values, L2 values]

% Collect the Label 1 files.
Files1 = {};
if MAPNfile
   Files1 = uipickfiles('FilterSpec', start_datadir, 'REFilter', ...
                        'MAPN.*\.mat', 'Prompt',                 ...
                        'Label 1 MAPN*.mat files');
else
   Files1 = uipickfiles('FilterSpec', start_datadir, 'REFilter', ...
                        '.*_Results.*\.mat', 'Prompt',           ...
                        'Label 1 _Results*.mat files');
end
n_files = numel(Files1);
if n_files == 0
   error('No files chosen for Label 1!');
end

% Collect the Label 2 files.
Files2 = {};
if MAPNfile
   Files2 = uipickfiles('FilterSpec', start_datadir, 'REFilter', ...
                        'MAPN.*\.mat', 'Prompt',                 ...
                        'Label 2 MAPN*.mat files');
else
   Files2 = uipickfiles('FilterSpec', start_datadir, 'REFilter', ...
                        '.*_Results.*\.mat', 'Prompt',           ...
                        'Label 2 _Results*.mat files');
end

if numel(Files2) ~= n_files
   error(['The two sets should each contain the same number of files!\n', ...
          'Label 1 (%d) versus Label 2 (%d)'], n_files, numel(Files2));
end

% Select the ROIs over all images.
n_ROIs_ALL = PA.defineROIs2(Files1, Files2, Pixel2nm, Color, ROI_sizes, ...
                            results_dir,                                ...
                            RegistrationNeeded, RegistrationViaDC);

%% ---------- Analyze the ROIS one-by-one

% NOTE:
%    results_c{1:n_files}{1:n_ROIs}{1:2}   clustering results by ROI and label
%    results_cs{1:n_files}{1:n_ROIs}       cluster separations between 2 labels
%    results_ls{1:n_files}{1:n_ROIs}       loc. separations between 2 labels
%    results_o1{1:n_files}{1:n_ROIs}       L1 clusters/L2 localizations overlap
%    results_o2{1:n_files}{1:n_ROIs}       L2 clusters/L1 localizations overlap
%    results_pcc{1:n_files}{1:n_ROIs}      pair cross-correlation results
%    resultsRC_pcc{1:n_files}              as above for ROIs combined
% The {1:n_files} is added at this level on the call to doAnalysis.

% Analyze the ROIs in each file.
opts = options;
if ~any(contains(opts, "plotting"))
   opts(end + 1) = "plotting";
end
files = dir(fullfile(data_dir, '*_ROIs.mat'));
n_files = numel(files);
fprintf('n_files = %d\n', n_files);
if n_files == 0
   error('No files chosen!');
end
%[~, files] = ...
%   smi_helpers.selectFiles(results_dir, '_ROIs.mat files', '*_ROIs.mat');
results_c     = cell(n_files, 1);
results_cs    = cell(n_files, 1);
results_ls    = cell(n_files, 1);
results_o1    = cell(n_files, 1);
results_o2    = cell(n_files, 1);
results_pcc   = cell(n_files, 1);
resultsRC_pcc = cell(n_files, 1);
for i = 1 : n_files
   filename = files(i).name;
   %filename = files{i};
   clear n_ROIs RoI
   % Load ROIs from file i.
   load(fullfile(data_dir, filename));
   % Clean up the description.
   desc = regexprep(filename, '_ROIs.mat', '');
   desc = regexprep(desc, '_ResultsStruct_BaGoL', '');
   desc = regexprep(desc, '_ResultsStruct', '');
   [results_pcc{i}, resultsRC_pcc{i}, results_c{i}, results_cs{i}, ...
    results_ls{i}, results_o1{i}, results_o2{i}] = ...
      PA.doAnalysis(n_ROIs, RoI, ROI_sizes, desc, particles, results_dir,  ...
                    opts, Pixel2nm, HistBinSize, RmaxAxisLimit, E, minPts, ...
                    PlotNonOverlap, Color);
end
analysis = 'analysis';
save(fullfile(results_dir, analysis), 'results_pcc', 'resultsRC_pcc', ...
     'results_c', 'results_cs', 'results_ls', 'results_o1', 'results_o2');
% Save analysis.mat
fprintf('Done one-by-one.\n');

%% ---------- Analyze a group of files of ROIs all together

%files = dir(fullfile(data_dir, '*_ROIs.mat'));
[~, files] = ...
   smi_helpers.selectFiles(data_dir, '_ROIs.mat files', '*_ROIs.mat');
n_files = numel(files);
fprintf('n_files = %d\n', n_files);
if n_files == 0
   error('No files chosen!');
end
answer = inputdlg('Name of combined results file:');
desc = answer{1};
if ~startsWith(desc, 'ALL')
   desc = ['ALL_', desc];
end

n_ROIs_ALL = 0;
for i = 1 : n_files
   %filename = files(i).name;
   filename = files{i};
   clear n_ROIs RoI
   % Load ROIs from file i.  Check for exclusions.
   load(fullfile(data_dir, filename));
   j = find(i == excludedFiles);
   if ~isempty(ExcludeROI) && ~isempty(j);
      excludedROIs = ExcludeROI{j}{2};
      RoI(excludedROIs) = [];
      n_ROIs = numel(RoI);
   end
   n_ROIs_ALL = n_ROIs_ALL + n_ROIs;
end
fprintf('n_ROIs_ALL = %d\n', n_ROIs_ALL);

% Combine all the ROIs together and reanalyze.
opts = options;
if ~any(contains(opts, "combined"))
   opts(end + 1) = "combined";
end
RoI_ALL = cell(n_ROIs_ALL, 1);
k = 0;
for i = 1 : n_files
   %filename = files(i).name;
   filename = files{i};
   clear n_ROIs RoI
   % Load ROIs from file i.  Check for exclusions.
   load(fullfile(data_dir, filename));
   j = find(i == excludedFiles);
   if ~isempty(ExcludeROI) && ~isempty(j);
      excludedROIs = ExcludeROI{j}{2};
      RoI(excludedROIs) = [];
      n_ROIs = numel(RoI);
   end
   for j = 1 : n_ROIs
      k = k + 1;
      RoI_ALL{k} = RoI{j};
   end
end

% Save analysisALL.mat
analysis = ['analysis', desc];
[results_pcc_ALL, resultsRC_pcc_ALL, results_c_ALL, results_cs_ALL, ...
 results_ls_ALL, results_o1_ALL, results_o2_ALL] = ...
   PA.doAnalysis(n_ROIs_ALL, RoI_ALL, ROI_sizes, desc, particles,         ...
                 results_dir, opts, Pixel2nm, HistBinSize, RmaxAxisLimit, ...
                 E, minPts, PlotNonOverlap, Color);
save(fullfile(results_dir, analysis), 'results_pcc_ALL',     ...
     'resultsRC_pcc_ALL', 'results_c_ALL', 'results_cs_ALL', ...
     'results_ls_ALL', 'results_o1_ALL', 'results_o2_ALL');
fprintf('Done %s.\n', desc);
