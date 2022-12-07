function hierBaGoL_run(Files, DataROI, Results_BaGoL, BaGoLParams, ROIs)
%hierBaGoL_run runs one or more BaGoL analyses.
% It is called by hierBaGoL_wrapper and calls hierBaGoL_analysis on each
% individual dataset to be processed, so acts as a dispatch intermediary.
% A single dataset is run directly, while a set of datasets are run in
% parallel using a parfor loop.
%
% INPUTS:
%    Files            cell array of full filepaths of SMLM datasets (containing
%                     one SMD structure per file) to be processed by BaGoL
%    DataROI          array of ROIs, one per file (pixels):
%                     [Xmin, Xmax, Ymin, Ymax]
%    Results_BaGoL    output directory name for BaGoL results; the results will
%                     be saved in SaveDir = fullfile(DataDir, Results_BaGoL),
%                     where DataDir is deduced from the input file path
%    BaGoLParams      see hierBaGoL_analysis for details
%    ROIs             for a single file of type _ROIs.mat, the ROIs selected
%                     via ClusterInterface.defineROIs using
%                     smi_helpers.ROITools.getROI are obtained
% OUTPUTS:
%    See hierBaGoL_analysis for details.

if exist('ROIs', 'var')
   ROIs = true;   % *_ROIs.mat file was input
else
   ROIs = false;
end

if ~iscell(Files)
   Files = { Files };
end
n_files = numel(Files);
fprintf('Files to analyze = %d\n', n_files);

if ~isempty(DataROI) && size(DataROI, 1) ~= n_files
   error('DataROI must either be [] or contain %d rows!', n_files);
end
% SE_Adjust is the standard error inflation (pixel), provided either as a
% constant (for all datasets) or an array, one per dataset
if numel(BaGoLParams.SE_Adjust) ~= 1 && numel(BaGoLParams.SE_Adjust) ~= n_files
   error('DataROI must contain 1 or %d values!', n_files);
end

if n_files > 0
   status = zeros(n_files, 1);

   if n_files == 1
      fprintf('(%d) %s ...\n', 1, Files{1});
      [DataDir, File, Ext] = fileparts(Files{1});
      SaveDir = fullfile(DataDir, Results_BaGoL);

      % Run hierBaGoL_analysis.
      %try
         data = load(Files{1});
         BaGoLParams.DataROI = DataROI;
         warning('OFF', 'stats:kmeans:FailedToConvergeRep');
         smi.BaGoL.hierBaGoL_analysis(data.SMD, File, SaveDir, BaGoLParams);
         status(1) = 1;
      %catch ME
      %   fprintf('### PROBLEM with %s ###\n', Files{i});
      %   fprintf('%s\n', ME.identifier);
      %   fprintf('%s\n', ME.message);
      %   status(i) = -1;
      %end
      fprintf('DONE (%d) %s.\n', 1, Files{1});
      fprintf('\n');
   else
      % When using parfor, dispatch the dataset analyses to as many workers as
      % available.
      delete(gcp('nocreate'));
      MachineInfo = parcluster();
      NumWorkers = MachineInfo.NumWorkers;
      parpool('local', min(NumWorkers, n_files));

      parfor i = 1 : n_files
         fprintf('(%d) %s ...\n', i, Files{i});
         [DataDir, File, Ext] = fileparts(Files{i});
         SaveDir = fullfile(DataDir, Results_BaGoL);

         % Set up BGLParams for parallel processing via a parfor loop.
         BGLParams = BaGoLParams;

         % Run hierBaGoL_analysis.
         try
            if ROIs
               filename = regexprep(Files{i}, '_ROI_[0-9][0-9]\.', '.');
               filename = regexprep(filename, ['Analysis', filesep], '');
               data = load(filename);
            else
               data = load(Files{i});
            end
            % If DataROI is defined, override the default value given in
            % BaGoLParams.
            if ~isempty(DataROI)
               BGLParams.DataROI = DataROI(i, :);
               fprintf( ...
                  'DataROI: [Xmin, Xmax, Ymin, Ymax] = [%g, %g, %g, %g]\n', ...
                       BGLParams.DataROI);
            end

            if numel(BaGoLParams.SE_Adjust) == 1          
               BGLParams.SE_Adjust = BaGoLParams.SE_Adjust;
            else
               BGLParams.SE_Adjust = BaGoLParams.SE_Adjust(i);
            end
            fprintf('SE_Adjust = %g\n', BGLParams.SE_Adjust);

            warning('OFF', 'stats:kmeans:FailedToConvergeRep');
            smi.BaGoL.hierBaGoL_analysis(data.SMD, File, SaveDir, BGLParams);
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
   end

   fprintf('BaGoL status by file (-1 PROBLEM, 0 NOT DONE, 1 DONE):\n');
   for i = 1 : n_files
      fprintf('[%2d] %2d   %s\n', i, status(i), Files{i});
   end
end
warning('ON', 'stats:kmeans:FailedToConvergeRep');
