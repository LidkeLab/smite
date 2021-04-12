% Script to produce BaGoL results from *_ResultsStruct.mat files.  The BaGoL
% results are placed in the subdirectory Results_BaGoL under the directory
% containing the *_ResultsStruct.mat files.

% Generic parameters
%start_DataDir = 'Y:\MJW\SR\play';
%start_DataDir = 'Y:\Sandeep\Genmab';
start_DataDir = 'Y:\Will K\Good Low Density IgE-AF647_dSTORM Data\Cell_02_b_test file\Label_01\No FC_Results';
%start_DataDir = '.';
BaGoLParams.ImageSize = 256;        % (pixel)
%BaGoLParams.PixelSize = 108.018;    % (nm) [TIRF]
BaGoLParams.PixelSize = 97.8;       % (nm) [sequential]
BaGoLParams.OutputPixelSize = 4;    %2; % pixel size for posterior images (nm)
BaGoLParams.N_Burnin1 = 3000;       % Length of Burn-in chain      (1st pass)
BaGoLParams.N_Trials1 = 3000;       % Length of post-burn-in chain (1st pass)
BaGoLParams.N_Burnin2 = 8000;       % Length of Burn-in chain      (2nd pass)
BaGoLParams.N_Trials2 = 4000;       % Length of post-burn-in chain (2nd pass)
% If Skip1stPass is true, skip the 1st pass in BaGoL_analysis and directly use
% the values for Lambda defined below.
BaGoLParams.Skip1stPass = true;
Results_BaGoL = 'Results_BaGoL';
% If Make2ndPassExpPrior is true, then the posterior Lambda from the 1st pass
% is converted into an exponential prior for the 2nd pass via:
%    [k2_prior, theta2_prior] = [1, k1_posterior * theta1_posterior]
% noting that k is the shape parameter (1 produces an exponential) and theta is
% the scale parameter for a gamma distribution, where the product k * theta is
% the mean of the distribution for the number of localizations per emitter, so
% the idea is to perserve this mean value while retaining an exponential shape.
% This is important in dSTORM.
BaGoLParams.Make2ndPassExpPrior = true;

% For very sparse datasets, turn off filtering (NNR = inf, NN = 0,
% NNNmax = inf).
%
% SE_Adjust adds to X_SE and Y_SE, so inflates the precision.  For DNA_PAINT
% data, SE_Adjust = 1--2 nm, while for dSTORM, slightly bigger values should be
% used.
%
% k and theta below are the shape and scale parameters for the Gamma
% probability distribution function.
%
% DataROI is defined when running BaGoL's first pass over only part of the
% image.  If DataROI is empty, use the whole image.
% 
BaGoLParams.IntensityCutoff = 5000; % Intensity cutoff
BaGoLParams.NNR = inf;              % At least NNN localizations within NNR(nm)
BaGoLParams.NNN = 0;                % Number of Nearest Neighbors
BaGoLParams.NNNmax = inf;           % No more than NNN locs within NNR (nm)
BaGoLParams.SE_Adjust = 3;          % Precision inflation applied to SE (nm)
BaGoLParams.ClusterDrift = 0;       % Expected magnitude of drift (nm/frame)
BaGoLParams.ROIsz = 500;            % ROI size for RJMCMC (nm)
%BaGoLParams.ROIsz = 200;            % ROI size for RJMCMC (nm)
BaGoLParams.OverLap = 50;           % Size of overlapping region (nm)
BaGoLParams.PreCluster = 50;        % Pre-clustering parameter (nm)
BaGoLParams.Lambda = [1, 1];        % [k, theta] parameters for gamma prior
BaGoLParams.DataROI = [];           % 1st pass [Xmin, Xmax, Ymin, Ymax] (pixel)
%BaGoLParams.DataROI = [100, 200, 100, 200];

% Sparse data parameters
%BaGoLParams.IntensityCutoff = 5000; % Intensity cutoff
%BaGoLParams.NNR = 8;                % At least NNN localizations within NNR(nm)
%BaGoLParams.NNN = 1;                % Number of Nearest Neighbors
%BaGoLParams.NNNmax = inf;           % No more than NNN locs within NNR (nm)
%BaGoLParams.SE_Adjust = 1;          % Precision inflation applied to SE (nm)
%BaGoLParams.ClusterDrift = 0;       % Expected magnitude of drift (nm/frame)
%BaGoLParams.ROIsz = 500;            % ROI size for RJMCMC (nm)
%BaGoLParams.OverLap = 50;           % Size of overlapping region (nm)
%BaGoLParams.PreCluster = 50;        % Pre-clustering parameter (nm)
%BaGoLParams.Lambda = [1.8, 4];      % [k, theta] parameters for gamma prior
%BaGoLParams.DataROI = [];           % 1st pass [Xmin, Xmax, Ymin, Ymax] (pixel)

% Dense data parameters
%BaGoLParams.IntensityCutoff = 5000; % Intensity cutoff
%BaGoLParams.NNR = 30;               % At least NNN localizations within NNR(nm)
%BaGoLParams.NNN = 1;                % Number of Nearest Neighbors
%BaGoLParams.NNNmax = inf;           % No more than NNN locs within NNR (nm)
%BaGoLParams.SE_Adjust = 1;          % Precision inflation applied to SE (nm)
%BaGoLParams.ClusterDrift = 0;       % Expected magnitude of drift (nm/frame)
%BaGoLParams.ROIsz = 1000;           % ROI size for RJMCMC (nm)
%BaGoLParams.OverLap = 50;           % Size of overlapping region (nm)
%BaGoLParams.PreCluster = 80;        % Pre-clustering parameter (nm)
%BaGoLParams.Lambda = [0.25, 200];   % [k, theta] parameters for gamma prior
%BaGoLParams.DataROI = [];           % 1st pass [Xmin, Xmax, Ymin, Ymax] (pixel)

c = SMA_Cluster();
Files = c.uipickfiles('FilterSpec', start_DataDir, 'REFilter', ...
                      '.*_ResultsStruct.*\.mat', 'Prompt',     ...
                      '_ResultsStruct.mat files');

if iscell(Files)
   n_files = numel(Files);
else
   n_files = 0;
end
fprintf('Files to analyze = %d\n', n_files);
if n_files > 0
   status = zeros(n_files, 1);

   delete(gcp('nocreate'));
   MachineInfo = parcluster();
   NumWorkers = MachineInfo.NumWorkers;
   parpool('local', min(NumWorkers, n_files));

   parfor i = 1 : n_files
      fprintf('(%d) %s ...\n', i, Files{i});
      [DataDir, File, Ext] = fileparts(Files{i});
      SaveDir = fullfile(DataDir, Results_BaGoL);
      try
         warning('OFF', 'stats:kmeans:FailedToConvergeRep');
         BaGoL_analysis(File, DataDir, SaveDir, BaGoLParams);
         status(i) = 1;
      catch ME
         fprintf('### PROBLEM with %s ###\n', Files{i});
         fprintf('%s\n', ME.identifier);
         fprintf('%s\n', ME.message);
         status(i) = -1;
      end
      fprintf('DONE (%d) %s.\n', i, Files{i});
      fprintf('\n');
   end

   fprintf('BaGoL status by file (-1 PROBLEM, 0 NOT DONE, 1 DONE):\n');
   for i = 1 : n_files
      fprintf('[%2d] %2d   %s\n', i, status(i), Files{i});
   end
end
warning('ON', 'stats:kmeans:FailedToConvergeRep');
fprintf('Done BaGoL.\n');
