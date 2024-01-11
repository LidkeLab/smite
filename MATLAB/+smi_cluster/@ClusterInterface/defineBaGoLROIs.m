function defineBaGoLROIs(pathnameR, filesR, pathnameB, filesB, MAPNfile)
% ---------- Possibly, define BaGoL ROIs from previous ROIs and BaGoL results
% Often, ROIs are defined from SR data, then need to be transferred to BaGoL
% (MAPN) processing of that SR data.  The BaGoL coordinates will replace the
% SR coordinates in the original ROI files.  Note that the two sets of files
% provided to this function should be in corresponding order.
%
% INPUTS:
%    pathnameR   path to where the ROI files below are located
%    filesR      _Results_ROIs.mat files to collect ROI info from
%    pathnameB   path to where the BaGoL files below are located
%    filesB      BaGoL files to collect ROI coordinate info from
%    MAPNfile    if false, assume BaGoL_Results_*_Results*.mat files,
%                otherwise if true, assume MAPN_*.mat files
%
% OUTPUTS:
%    Saves pathnameB/Analysis/*_BaGoL_ROIs.mat

% Created by
%    Michael J. Wester (2022)

   n_files = numel(filesR);
   if n_files ~= numel(filesB)
      error('n_filesROIs != n_filesBaGoL!')
   end
   results_dir = fullfile(pathnameB, 'Analysis');
   % Create results_dir if it does not already exist.
   if ~isfolder(results_dir)
      mkdir(results_dir);
   end

   for l = 1 : n_files
      dataR = load(fullfile(pathnameR, filesR{l}));
      dataB = load(fullfile(pathnameB, filesB{l}));
      if MAPNfile
         dataB_MAPN = dataB.MAPN;
      else
         dataB_MAPN = dataB.BGL.MAPN;
      end
      %[~, fileB, ~] = fileparts(dataR.ResultsFile);
      % dataR.ResultsFile cannot be relied on for old runs, so this is safer.
      [~, fileB, ~] = fileparts(filesR{l});
      short = regexprep(fileB, '.mat$', '');
      ResultsFile = dataR.ResultsFile;
      Pixel2nm = dataR.Pixel2nm;
      XYsize = dataR.XYsize;
      n_ROIs = dataR.n_ROIs;
      RoI = cell(1, n_ROIs);
      for i = 1 : n_ROIs
         ROI = dataR.RoI{i}.ROI;
         RoI{i}.ROI = ROI;
         xmin = ROI(1);   xmax = ROI(2);   ymin = ROI(3);   ymax = ROI(4);
         n_labels = numel(dataR.RoI{1}.X);
         for j = 1 : n_labels
            Xnm = dataB_MAPN.X; % * Pixel2nmGlobal;
%           if RT.GaussIm    % if GaussIm, have DIPimage style coordinates
%              Ynm = XYsize(2) - dataB_MAPN.Y;
%           else
               Ynm = dataB_MAPN.Y;
%           end
            k = xmin <= Xnm & Xnm <= xmax & ymin <= Ynm & Ynm <= ymax;
            fprintf('File %d ROI %d Label %d: %d points\n', l, i, j, sum(k));
            % Check that k includes some points!
            if any(k)
               % Sanity check!
               if min(Ynm(k)) < ymin || max(Ynm(k)) > ymax
                  error('ROI -> BaGoL ROI: ROI = %f %f, min/max Ynm = %f %f',...
                        ROI(3), ROI(4), min(Ynm(k)), max(Ynm(k)));
               end

               RoI{i}.X{j}     = Xnm(k);
               RoI{i}.Y{j}     = Ynm(k);
               RoI{i}.X_STD{j} = dataB_MAPN.X_SE(k);
               RoI{i}.Y_STD{j} = dataB_MAPN.Y_SE(k);
            else
               warning('ROI -> BaGoL ROI: no points!');
               RoI{i}.X{j}     = [];
               RoI{i}.Y{j}     = [];
               RoI{i}.X_STD{j} = [];
               RoI{i}.Y_STD{j} = [];
            end
         end
      end
      save(fullfile(results_dir, [short, '_BaGoL_ROIs.mat']), ...
           'ResultsFile', 'Pixel2nm', 'XYsize', 'n_ROIs', 'RoI');
   end
   if MAPNfile
      txt = 'MAPN ROIs';
   else
      txt = 'ROIs';
   end
   fprintf('Done BaGoL %s: %s.\n', txt, results_dir);

end
