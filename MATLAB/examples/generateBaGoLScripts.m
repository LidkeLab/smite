% Generate a series of self-contained scripts to run BaGoL over a number of
% cells, one script per cell.  This allows the BaGoL analysis of each cell to
% be run independently, which is useful when running on different servers
% sharing a common data space.  The main parameters of this script generation
% process are
% -  the tenplate location (and the template itself, of course---see below),
% -  the directory of SMITE results, which is assumed to to have a subdirectory
%    ResultsStructs containing individual "Cell_*_Results.mat" files.  The cell
%    number is embedded in the name of the _Results.mat files,
% -  a short name for the BaGoL output directory (Results_BaGoL) to which the
%    cell number will be appended, e,g., B12_20240816_02,
% -  the fullpath of the BaGoL directory that will contain the scripts
%    produced.
clear all

% BaGoL template script location.  Note that 3 lines in this file start with a
% +.  These lines will be replaced by the appropriate information
% (BaGoL_Results directory, D1 data directory, and cell data filename) when
% this script is run.
Template = '/mnt/nas/cellpath/Genmab/Data/Analyses/BaGoL_TEMPLATE.m';
ResultsDir = '/mnt/nas/cellpath/Genmab/Data/2024-05-17_HeLa_SaturatingIgG10min/HeLa_IgG1-B12-AF647_5ugml_10min/Results';
Results_BaGoL = 'B12_20240816';
BaGoLDir = fullfile(ResultsDir, 'BGL');

cells = dir(fullfile(ResultsDir, 'ResultsStructs', 'Cell_*_Results.mat'));
n_cells = numel(cells);
if n_cells == 0
   error('No cells found!');
end

if ~exist(BaGoLDir, 'dir')
   mkdir(BaGoLDir);
end

% For each Cell_*_Results.mat file found in ResultsStructs ...
for i = 1 : n_cells
   cellname   = cells(i).name;
   cellnumber = regexprep(cellname, '^Cell_([0-9][0-9]).*mat$', '$1');
   out_name   = fullfile(BaGoLDir, sprintf('Cell%s', cellnumber));

   % Copy the template to the script, making substitutions on lines beginning
   % with a +.
   in  = fopen(Template, 'r');
   out = fopen(out_name, 'w');

   while true
      textline = fgetl(in);
      if isempty(textline)
         fprintf(out, '\n');
      elseif textline(1) == -1   % EOF indication
         break;
      elseif textline(1) == '+'
         % BaGoL output directory name.
         if startsWith(textline, '+Results_BaGoL')
            fprintf(out, 'Results_BaGoL = ''%s_%s'';\n', ...
                         Results_BaGoL, cellnumber);
         % SMITE Results directory location.
         elseif startsWith(textline, '+D1')
            fprintf(out, 'D1 = ''%s'';\n', ...
                         fullfile(ResultsDir, 'ResultsStructs'));
         % SMITE Cell_*_Results file.
         elseif startsWith(textline, '+fullfile(D1,')
            fprintf(out, 'fullfile(D1, ''%s'');\n', cellname);
         else
         end
      else
         fprintf(out, '%s\n', textline);
      end
   end

   fclose(in);
   fclose(out);
end
