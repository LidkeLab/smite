function combineBaGoLROIs(pathnameR, filesR, pathnameB, filesB, MAPNfile, ...
                          keep_numbering, GaussIm)
% ---------- Possibly, combine individually processed BaGoL ROIs into a single
%            _ROIs.mat file using the _ROIs.mat file that was used to define
%            the ROIs originally from the SR data
% If datasets are too dense, it becomes necessary to process each ROI
% separately.  This function will combine multiple BaGoL processed ROIs from
% multiple (biological) cells into appropriately named combined _ROIs.mat
% files.  This assumes the various files involved are following a naming
% convention.  For example, the original collection of ROIs produced from the
% SR data (combining both labels):
%    Cell_01_Label_01_Results_ROIs.mat
%    Cell_02_Label_01_Results_ROIs.mat
%    ...
% The BaGoL MAPN coordinates per ROI are named:
%    MAPN_Cell_01_Label_01_Results_ROI_01.mat
%    MAPN_Cell_01_Label_01_Results_ROI_02.mat
%    ...
%    MAPN_Cell_01_Label_02_Results_ROI_01.mat
%    ...
%    MAPN_Cell_02_Label_01_Results_ROI_01.mat
%    ...
% The collection of ROIs (combining both labels) produced from the BaGoL data:
%    Cell_01_Results_BaGoL_ROIs.mat
%    Cell_02_Results_BaGoL_ROIs.mat
%    ...
%
% INPUTS:
%    pathnameR   path to where the ROIs file below are located
%    filesR      _Results_ROIs.mat files to collect ROI info from
%    pathnameB   path to where the BaGoL files below are located
%    filesB      BaGoL files to collect ROI coordinate info from
%    MAPNfile    if false, assume BaGoL_Results_*_Results*.mat files,
%                otherwise if true, assume MAPN_*.mat files
%    keep_numbering retain the ROI numbering even if there are missing ROIs
%                (which will be treated as empty) [default: false]
%    GaussIm     true if the original SR ROIs were chosen from a Gaussian image
%                [default: false]
%
% OUTPUTS:
%    Saves pathnameB/Analysis/*_BaGoL_ROIs.mat

% Created by
%    Michael J. Wester (2022)

   if ~exist('keep_numbering', 'var')
      keep_numbering = false;
   end
   if ~exist('GaussIm', 'var')
      GaussIm = false;
   end

   results_dir = pathnameB;

   n_files = numel(filesB);
   n_cells = numel(filesR);
   k = 0;
   for cc = 1 : n_cells
      % Extract the cell number from the file name.
      c = str2num(regexprep(filesR{cc}, '^.*Cell_([0-9][0-9]).*$', '$1'));
      cellno(c) = cc;
      % Copy the global parameters from the original SR _ROIs.mat file.
      dataR = load(fullfile(pathnameR, filesR{cc}));
      if isfield(dataR, 'ResultsFile');
         ResultsFile = dataR.ResultsFile;
      else
         ResultsFile = dataR.Label1;
         ResultsFile = regexprep(ResultsFile, '_Label_[0-9][0-9]', '');
         ResultsFile = regexprep(ResultsFile, '_L[0-9][0-9]', '');
      end
      if isfield(dataR, 'Label1');
         Label1 = dataR.Label1;
      end
      if isfield(dataR, 'Label2');
         Label2 = dataR.Label2;
      end
      Pixel2nm = dataR.Pixel2nm;
      XYsize = dataR.XYsize;
      n_labels = numel(dataR.RoI{1}.X);
      n_ROIs = dataR.n_ROIs;

      % Consistency check.
      %n_files = numel(filesB);
      %if n_labels * n_ROIs ~= numel(filesB)
      %   error('n_labels * n_ROIs != n_filesBaGoL!')
      %end

      % Gather the coordinates and standard errors from the individual
      %_ROI_*.mat files containing BaGoL processed data.
      found = [];   % a list of found ROIs
      RoI = cell(n_ROIs, 1);
      for jj = 1 : n_labels
         for ii = 1 : n_ROIs
            k = k + 1;
            if k <= n_files
               fileB = filesB{k};
               dataB = load(fullfile(pathnameB, fileB));
               % Extract the ROI number and label number from the file name. 
               i = str2num(regexprep(fileB, '^.*ROI_([0-9][0-9]).*$',   '$1'));
               j = str2num(regexprep(fileB, '^.*Label_([0-9][0-9]).*$', '$1'));
               d = str2num(regexprep(fileB, '^.*Cell_([0-9][0-9]).*$',  '$1'));
               if d ~= c
                  error('Inconsistency in Cell #: expected %d, got %d', c, d);
               end
               found = [found, i]; 
               
               RoI{i}.ROI = dataR.RoI{i}.ROI;
               if MAPNfile
                  dataB_MAPN = dataB.MAPN;
               else
                  dataB_MAPN = dataB.BGL.MAPN;
               end
               RoI{i}.X{j} = dataB_MAPN.X;
               RoI{i}.Y{j} = dataB_MAPN.Y;
               RoI{i}.X_STD{j} = dataB_MAPN.X_SE;
               RoI{i}.Y_STD{j} = dataB_MAPN.Y_SE;
            elseif keep_numbering
               % Not enough BaGoL files; insert empty data for the missing ROIs
               % to keep ROI numbering consistent.
               missing = setdiff(1 : n_ROIs, found);
               for kk = 1 : numel(missing)
                  RoI{missing(kk)}.ROI = dataR.RoI{missing(kk)}.ROI;
                  RoI{missing(kk)}.X = {[],[]};
                  RoI{missing(kk)}.Y = {[],[]};
                  RoI{missing(kk)}.X_SE = {[],[]};
                  RoI{missing(kk)}.Y_SE = {[],[]};
               end
               break
            end
         end
      end

      % Eliminate empty ROIs.
      missing = setdiff(1 : n_ROIs, found);
      if ~keep_numbering && isempty(missing)
         RoI_new = cell(n_ROIs - numel(missing), 1);
         k = 1;
         for ii = 1 : n_ROIs
            RoI_new{k} = RoI{ii};
            if ~ismember(ii, missing)
               k = k + 1;
            end
         end
         n_ROIs = n_ROIs - numel(missing); 
         RoI = RoI_new;
      end

      [~, fileB, ~] = fileparts(ResultsFile);
      short = regexprep(fileB, '.mat$', '');

      for i = 1 : n_ROIs
         ROI = RoI{i}.ROI;
         xmin = ROI(1);   xmax = ROI(2);   ymin = ROI(3);   ymax = ROI(4);
         for j = 1 : n_labels
            Xnm = RoI{i}.X{j};
%           if RT.GaussIm    % if GaussIm, have DIPimage style coordinates
%              Ynm = XYsize(2) - RoI{i}.Y{j};
%           else
               Ynm = RoI{i}.Y{j};
%           end
            % If the original SR ROIs were chosen via GaussIm and then run
            % individually through BaGoL, ymin and ymax are flipped about the
            % center of the image, so need to be corrected.  This problem
            % really should be fixed when the BaGoL ROI is defined, but for now
            % ...
            if GaussIm && (ymin > min(Ynm) || ymax < max(Ynm))
               ymin = XYsize(2) - ymax;
               ymax = XYsize(2) - ymin;
            end
            l = xmin <= Xnm & Xnm <= xmax & ymin <= Ynm & Ynm <= ymax;
            fprintf('Cell %d ROI %d Label %d: %d points\n', c, i, j, sum(l));
            % Check that l includes some points!
            if any(l)
               % Sanity check!
               if min(Ynm(l)) < ymin || max(Ynm(l)) > ymax
                  error('ROI -> BaGoL ROI: ROI = %f %f, min/max Ynm = %f %f',...
                        ROI(3), ROI(4), min(Ynm(l)), max(Ynm(l)));
               end
            else
               warning('ROI -> BaGoL ROI: no points!');
            end
         end
      end

      fprintf('Saving %s ...\n', [short, '_BaGoL_ROIs.mat']);
      if exist('ResultsFile', 'var')
         save(fullfile(results_dir, [short, '_BaGoL_ROIs.mat']), ...
              'ResultsFile', 'Pixel2nm', 'XYsize', 'n_ROIs', 'RoI');
      else
         save(fullfile(results_dir, [short, '_BaGoL_ROIs.mat']), ...
              'Label1', 'Label2', 'Pixel2nm', 'XYsize', 'n_ROIs', 'RoI');
      end
   end % for cc

   if MAPNfile
      txt = 'MAPN ROIs';
   else
      txt = 'ROIs';
   end
   fprintf('Done combining BaGoL %s: %s.\n', txt, results_dir);

end
