% Perform cluster analyses for comparison of experimental conditions.  This is
% a batch version of simpleROIcluster for combining statistics for one or more
% conditions. *_results.mat files produced by simpleROIcluster (Statistics for
% a single condition) or singleConditionDriver are the inputs to these
% analyses.
%
% NOTE: smite and PlotSpread should be on your MATLAB path.
%
% Work flow:
%   Set important parameters:
%      ROI_sizes, Pixel2nm, start_datadir, oneROI
%   Combined statistics for one or more conditions:
%            multiple [condition]_results_mat -> 
%      multiple conditions are plotted on the same graph and placed in a
%      subdirectory for this particular analysis.  The filenames are assumed
%      to have the structure:
%         experimentalConditions#analysisConditions_results.mat
%      where experimentalConditions contains no #.  Statistics for a single
%      condition (see simpleROIcluster) produces filenames in this format, so
%      can be copied directly.
%      A typical example of a *_results.mat file is:
%         ConditionA#DBSCAN_N=3,E=50_results.mat
%   ALTERNATIVELY, combined statistics for multiple conditions and experiments
%      (as above, where line colors and types for CDF2 plots are set by hand)
%      See simpleROIcluster (or in the future combinedStats2).
%
% Avoid spaces in filenames and conditions!!!

%% ---------- Set important parameters

doHopkins = true;   % Hopkins' test can be time consuming for dense ROIs
% Set to false for very dense ROIs to avoid crashes due to lack of memory.
doSigmaActual = false;

ROI_sizes = [3000, 3000];   % [delta_x, delta_y] (nm)
A_ROI = prod(ROI_sizes);    % ROI area (nm^2)
%Pixel2nm = 16000/150;       % conversion factor from pixels to nm
%Pixel2nm = 108.018;         % pixels to nm [TIRF]
Pixel2nm = 97.8;            % pixels to nm [sequential]
Pixel2nmGlobal = Pixel2nm;

% ClusterInterface contains various helper routines used by this script.
CI = smi_cluster.ClusterInterface();

% Often, for BaGoL analyses, it is simpler to use a single, large, encompassing
% ROI rather than a series of small ROIs.
oneROI = false;
if oneROI
   ROI_sizes = [256, 256] * Pixel2nm;   % (nm)
   A_ROI = prod(ROI_sizes);
end

% Select the files starting from start_datadir.
start_datadir = '.';

fprintf('Done set parameters.\n');

%% ---------- Combined statistics for multiple conditions and experiments

SC = smi_cluster.StatisticsClustering();
% Make various plots:
%    'f'   frequency
%    'n'   normalized
%    'p'   PDF
%    'c'   CDF
%    'C'   CDF (alternative)
%    's'   PlotSpread
%    'S'   PlotSpread (bars for mean & median)
%    'x'   box
%    'b'   bar
SC.PlotDo = 'CSx';
% Red mean, green median (2 only mean, 3 only median) for PlotSpread plots.
SC.ShowMM = 1;
% Options for CDF2 plots are: 'plot', 'semilogx', 'semilogy', 'loglog'.
SC.LinLog = 'semilogx';
SC.Ylim = [0.01, 1];

% Use the same color for each experiment under the same condition, different
% colors for different conditions.  There is one entry per file below.
colors = ['b', 'b', 'r', 'r', 'g', 'g'];

% Use different line types to distingush same colored lines when desired.
line_type = {'-', '--', '-', '--', '-', '--'};

%[pathname, files] = smi_helpers.selectFiles(start_datadir, ...
%                       '_results.mat files', '*_results.mat');
%answer = inputdlg('Output directory identifier:');
%base_name = answer{1};

% Manual specification: one *_results.mat file per condition and experiment,
% all collected in the directory pathname.
pathname = start_datadir;
N = 3;
E = 30;
files = {
sprintf('ConditionA_Exp1#DBSCAN_N=%d,E=%d_results.mat', N, E)
sprintf('ConditionA_Exp2#DBSCAN_N=%d,E=%d_results.mat', N, E)
sprintf('ConditionB_Exp1#DBSCAN_N=%d,E=%d_results.mat', N, E)
sprintf('ConditionB_Exp2#DBSCAN_N=%d,E=%d_results.mat', N, E)
sprintf('ConditionC_Exp1#DBSCAN_N=%d,E=%d_results.mat', N, E)
sprintf('ConditionC_Exp2#DBSCAN_N=%d,E=%d_results.mat', N, E)
};
base_name = 'CONSOLIDATED';

CI.combinedStatistics2(SC, colors, line_type, pathname, files, base_name, ...
                       A_ROI, doHopkins);
