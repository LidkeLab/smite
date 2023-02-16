function combineResults(Files, analysis_dir, out_file)
% ---------- Combine results from 1 or more conditions for further processing
% Multiple instantiations of a single condition potentially located in multiple
% directories are collected into a combined instance in a specified directory.
%
% INPUTS:
%    Files          _results.mat files to collect analyses from
%    analysis_dir   directory where to put the combined result
%    out_file       name of combined ouput _results.mat file
%
% OUTPUTS:
%    analysis_dir/out_file is saved where analysis_dir is typically
%    some_path/Analysis and outfile ends in _results.mat

% Created by
%    Michael J. Wester (2022)

   n_files = numel(Files);
   if n_files == 1
      copyfile(Files{1}, fullfile(analysis_dir, out_file));
   else
      if ~endsWith(out_file, '_results.mat')
         out_file = [out_file, '_results.mat'];
      end
      files = {};
      results = {};
      RoI = {};
      n_ROIs = 0;
      collected = {};
      for j = 1 : n_files
         cdata = load(Files{j});

         files = {files{:}, cdata.files{:}};
         results = {results{:}, cdata.results{:}};
         RoI = {RoI{:}, cdata.RoI{:}};

         % Pixel2nm should be consistent between result files.
         if j == 1
            Pixel2nm = cdata.Pixel2nm;
         elseif cdata.Pixel2nm ~= Pixel2nm
            error('Pixel2nm inconsistent: %g %g', Pixel2nm, cdata.Pixel2nm);
         end

         n_ROIs = n_ROIs + cdata.n_ROIs;

         fn = fieldnames(cdata.collected);
         for i = 1 : numel(fn)
            fn_i = fn{i};
            if j == 1
               if strcmp(fn_i, 'nC')
                  collected.nC = 0;
               else
                  collected.(fn_i) = [];
               end
            end
            cdatum = cdata.collected.(fn_i);
            if isscalar(cdatum)
               if strcmp(fn_i, 'nC')
                  collected.nC = collected.nC + cdatum;
               else
                  collected.(fn_i) = cdatum;
               end
            else
               collected.(fn_i) = [collected.(fn_i), cdatum];
            end
         end
      end
      save(fullfile(analysis_dir, out_file), 'files', 'Pixel2nm', 'n_ROIs', ...
                    'RoI', 'results', 'collected');
   end

   fprintf('Done.\n');

end
