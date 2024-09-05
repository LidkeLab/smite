% BaGoL run TEMPLATE file used by generateBaGolScripts.m

% Output directory name.
+Results_BaGoL = 'DF3_01';

% Generic parameters.
BaGoLParams.ImageSize = 256;        % (pixel)
%BaGoLParams.PixelSize = 108.018;    % (nm) [TIRF]
BaGoLParams.PixelSize = 97.8;       % (nm) [sequential]
BaGoLParams.OutputPixelSize = 4;    %2; % pixel size for posterior images (nm)
BaGoLParams.N_Burnin = 32000;       % Length of Burn-in chain
BaGoLParams.N_Trials = 8000;        % Length of post-burn-in chain
%BaGoLParams.N_Burnin = 64000;       % Length of Burn-in chain
%BaGoLParams.N_Trials = 16000;        % Length of post-burn-in chain
%BaGoLParams.N_Burnin = 8000;        % Length of Burn-in chain
%BaGoLParams.N_Trials = 4000;        % Length of post-burn-in chain
BaGoLParams.NSamples = 10;          % Number of samples before sampling Xi
BaGoLParams.ClusterDrift = 0;       % Expected magnitude of drift (nm/frame)

% Y_Adjust is sometimes needed to deal with lower left versus upper left
% y-origin issues.  Lower left with y increasing upwards is the default,
% requiring no changes, while upper left with y increasing downwards can
% sometimes occur, so Y must be replaced by Y - Y_Adjust, where Y_Adjust is the
% image size (see below) [pixels].
%BaGoLParams.Y_Adjust = BaGoLParams.ImageSize;
BaGoLParams.Y_Adjust = [];

% SE_Adjust adds to X_SE and Y_SE, so inflates the precision.  For DNA_PAINT
% data, SE_Adjust = 1--2 nm, while for dSTORM, slightly bigger values should
% be used.  Note that this quantity can be specified as an array of length
% n_files if applied differently to each file.
BaGoLParams.SE_Adjust = 5;          % Precision inflation applied to SE (nm)
%BaGoLParams.SE_Adjust = [0, 0];     % Precision inflation applied to SE (nm)

% The values for ROIsz and OverLap directly below are good for denser data as
% less computational effort is required, so the code runs faster.  The second
% set of values can be used for sparser data to generate larger ROIs, but may
% produce artifacts with dense data.
BaGoLParams.ROIsz = 500;            % ROI size for RJMCMC (nm)
BaGoLParams.OverLap = 50;           % Size of overlapping region (nm)
BaGoLParams.Cutoff = 25;            % Pre-clustering cutoff (nm)
%BaGoLParams.ROIsz = 50;             % ROI size for RJMCMC (nm)
%BaGoLParams.OverLap = 10;           % Size of overlapping region (nm)
%BaGoLParams.ROIsz = 100;            % ROI size for RJMCMC (nm)
%BaGoLParams.OverLap = 25;           % Size of overlapping region (nm)
%BaGoLParams.ROIsz = 500;            % ROI size for RJMCMC (nm)
%BaGoLParams.OverLap = 50;           % Size of overlapping region (nm)

% k and theta below are the shape and scale parameters for the Gamma
% probability distribution function.  If just one parameter is provided,
% a Poisson distribution is used.
BaGoLParams.Xi = [1, 10];           % [k, theta] parameters for gamma prior

% Note for batch runs, in which Files and DataROI are input by hand, please see
% ### comments below.
BaGoLParams.DataROI = [];           % [Xmin, Xmax, Ymin, Ymax] (pixel)
DataROI = [];

% If ROIs is true, the input file has ROIs already defined (*_ROIs.mat), so use
% them below if only one filename is provided.
%ROIs = true;
ROIs = false;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%start_DataDir = '.';
%Files = uipickfiles('FilterSpec', start_DataDir, 'REFilter', ...
%                    '.*_ResultsStruct.*\.mat', 'Prompt',     ...
%                    '_ResultsStruct.mat files');
% ### Comment out the 4 lines above and use the commented out lines below when
% ### making batch runs, for example, on CARC.  Here, the files to process are
% ### defined relative to the directory where hierBaGoL_wrapper is run.
% ### Absolute pathnames are also fine, especially when used in conjunction
% ### with fullfile.
+D1 = '/mnt/nas/cellpath/Personal Folders/Rachel/20231010_IgE-AntigenStrain/BaGoLHierCM/2024WillRedo/DF3/Results/ResultsStructs';
Files = {
+fullfile(D1, 'Cell_01_Label_01_Data_2021-3-26-23-24-49_Results.mat');
};

% DataROI is defined when running BaGoL over only part of the image.
% If DataROI is empty, use the whole image.
% 
% Define a single region of interest for each dataset (units are pixels).
% [YStart, XStart, YEnd, XEnd] = [163, 385, 233, 455]
% [Xmin, Xmax, Ymin, Ymax] (pixel)
% [385, 455, 163, 233]
%DataROI = [
%[120, 136, 190, 206]
%[110, 126,  90, 106]
%];

% Special case in which a single file is broken up into pre-selected ROIs if
% the variable ROIs is set to true.  The ROIs are processed via a parfor loop
% in smi.BaGoL.hierBaGoL_run.  See also comments at the top of this function.
if numel(Files) == 1 && ROIs
   [DataDir, File, Ext] = fileparts(Files{1});
   basename = fullfile(DataDir, 'Analysis', File);
   % Assume SMD files are of the form Cell_nn_Label_0n_Results.mat and ROI
   % files are of the form Cell_nn_Label_01_Results_ROIs.mat (and are located
   % in the subdirectory 'Analysis' of DataDir).
   filename = regexprep(basename, 'Label_02', 'Label_01');
 

   ROIsFile = load([filename, '_ROIs.mat']);
   n_ROIs = numel(ROIsFile.RoI);
   DataROI = zeros(n_ROIs, 4);
   for i = 1 : n_ROIs
      DataROI(i, :) = ROIsFile.RoI{i}.ROI ./ BaGoLParams.PixelSize;
   end

   Files = cell(n_ROIs, 1);
   for i = 1 : n_ROIs
      Files{i} = sprintf('%s_ROI_%02d.mat', basename, i);
   end

end

% ----------------------------------------------------------------------

% Run the BaGoL analyses.
smi.BaGoL.hierBaGoL_run(Files, DataROI, Results_BaGoL, BaGoLParams, ROIs);

fprintf('Done BaGoL.\n');
