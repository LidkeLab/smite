function combinedStatistics1(SC, pathname, files, base_name, A_ROI, doHopkins)
% ---------- Combined statistics for one or more conditions
% For multiple conditions, experimental (e.g., resting vs. activated) and/or
% analytical (e.g., various DBSCAN parameter combinations), produce plots
% that compare the results for a variety of studies (e.g., Hopkins' statistic,
% number of clusters per ROI, circular equivalent cluster radii, etc.).
%
% INPUTS:
%    SC          smi_cluster.StatisticsClustering object setting values for
%                properties PlotDo, ShowMM, LinLog, Ylim
%
%                SC = smi_cluster.StatisticsClustering();
%                % Make various plots:
%                %    'f'   frequency
%                %    'n'   normalized
%                %    'p'   PDF
%                %    'c'   CDF
%                %    'C'   CDF (alternative)
%                %    's'   PlotSpread
%                %    'S'   PlotSpread (bars for mean & median)
%                %    'x'   box
%                %    'b'   bar
%                SC.PlotDo = 'CSx';
%                % Red mean, green median (2 only mean, 3 only median)
%                % for PlotSpread plots.
%                SC.ShowMM = 1;
%                % Options for CDF2 plots are:
%                %    'plot', 'semilogx', 'semilogy', 'loglog'.
%                SC.LinLog = 'semilogx';
%                SC.Ylim = [0.01, 1];
%                % Produce a comma separated value file with the data for easy
%                % import into other software.
%                SC.CSV = true;
%
%                NOTE: SC.CSV is turned on and off internally to produce comma
%                separated values files for some plots.
%
%    pathname    path to where the single condition results are located
%    files       _results.mat files to collect info from
%    base_name   analysis identification string
%    A_ROI       ROI area (nm^2)
%    doHopkins   [DEFAULT = true] Hopkins' test can be time consuming for dense
%                ROIs
%
% OUTPUTS:
%    A variety of plots are produced via SC.plotCombined and are placed in a
%    subdirectory of pathname.  The subdirectory's name is derived from the
%    the analysisCondition taken from the first input filename (all filenames
%    compared are assumed to have been analyzed similarly), which is assumed
%    to have the form:
%       experimentalConditions#analysisConditions_results.mat
%    where experimentalConditions contains no #.  This function will extract
%    the various experimentalConditions being compared and use them to label
%    the plots created.

% Created by
%    Michael J. Wester (2022)

CSVglobal = false;

c = smi_cluster.Clustering();

n_files = numel(files);

% Make a subdirectory for this particular analysis.  The filenames are assumed
% to have the structure:
%    experimentalConditions#analysisConditions_results.mat
% where experimentalConditions contains no #.
econd = cell(1, n_files);
for i = 1 : n_files
   econd{i} = regexprep(files{i}, '^(.*)#.*_results.mat', '$1');
end
analysis = regexprep(files{1}, '^.*#(.*)_results.mat', '$1');
Cresults_dir = fullfile(pathname, [base_name, '_', analysis]);
SC.ResultsDir = Cresults_dir;
% Create Cresults_dir if it does not already exist.
if ~isfolder(Cresults_dir)
   mkdir(Cresults_dir);
end

SC.BaseName = '';

for j = 1 : n_files
   rdata{j} = load(fullfile(pathname, files{j}));
end

% Create file and format for P-values.
out = fopen(fullfile(Cresults_dir, 'P_values.txt'), 'w');
fmt = '';
for j = 1 : n_files
   fmt = [fmt, ' %.5e'];
end
fmt = [fmt, '\n'];

% Filter out extrema.
for j = 1 : n_files
%  rdata{j}.collected.n_clusts      = ...
%     rdata{j}.collected.n_clusts(rdata{j}.collected.n_clusts         < 400);
%  rdata{j}.collected.c_density     = ...
%     rdata{j}.collected.c_density(rdata{j}.collected.c_density       < 1.6e-5);
%  rdata{j}.collected.n_points      = ...
%     rdata{j}.collected.n_points(rdata{j}.collected.n_points         < 10000);
%  rdata{j}.collected.n_pts         = ...
%     rdata{j}.collected.n_pts(rdata{j}.collected.n_pts               < 600);
%  rdata{j}.collected.nn_dists      = ...
%     rdata{j}.collected.nn_dists(rdata{j}.collected.nn_dists         < 800);
%  rdata{j}.collected.equiv_radii   = ...
%     rdata{j}.collected.equiv_radii(rdata{j}.collected.equiv_radii   < 300);
%  rdata{j}.collected.nn_loc_dists  = ...
%     rdata{j}.collected.nn_loc_dists(rdata{j}.collected.nn_loc_dists < 400);
end

% Hopkins' statistic per ROI.
if doHopkins
   %SC.CSV = true;
   P_hopkins_ROI = ...
   SC.plotCombined(arrayfun(@(i) rdata{i}.collected.H_ROI, 1 : n_files, ...
                            'UniformOutput', false),                    ...
                   0.01, 'Hopkins per ROI', econd, 'hopkins_ROI');
   %SC.CSV = false;
   fprintf(out, 'P_hopkins per ROI =\n');
   fprintf(out, fmt, P_hopkins_ROI);
   fprintf(out, '\n');
end

% Number of clusters per ROI.
P_n_clusts = ...
SC.plotCombined(arrayfun(@(i) rdata{i}.collected.n_clusts, 1 : n_files,  ...
                         'UniformOutput', false),                        ...
                1, '# of clusters per ROI', econd, 'n_clusts_ROI');
fprintf(out, 'P_# of clusters per ROI =\n');
fprintf(out, fmt, P_n_clusts);
fprintf(out, '\n');

% Density of clusters per ROI.
P_c_density = ...
SC.plotCombined(arrayfun(@(i) rdata{i}.collected.c_density*1000^2,  ...
                              1 : n_files, 'UniformOutput', false), ...
                1/(A_ROI/1000^2),                                   ...
                'density of clusters per ROI (\mu m^{-2})',         ...
                econd, 'c_density_ROI');
fprintf(out, 'P_density of clusters per ROI =\n');
fprintf(out, fmt, P_c_density);
fprintf(out, '\n');

% Number of localizations per ROI.
if CSVglobal
   SC.CSV = true;
end
P_n_local = ...
SC.plotCombined(arrayfun(@(i) rdata{i}.collected.n_points, 1 : n_files, ...
                         'UniformOutput', false),                       ...
                [], '# of localizations per ROI', econd, 'n_local_ROI');
SC.CSV = false;
fprintf(out, 'P_# of localizations per ROI =\n');
fprintf(out, fmt, P_n_local);
fprintf(out, '\n');

% Density of localizations per ROI.      
P_den_local = ...
SC.plotCombined(arrayfun(@(i) rdata{i}.collected.n_points/A_ROI*1000^2, ...
                              1 : n_files, 'UniformOutput', false),     ...
                [], 'density of localizations per ROI (\mu m^{-2})',    ...
                econd, 'den_local_ROI');
fprintf(out, 'P_density of localizations per ROI =\n');
fprintf(out, fmt, P_den_local);
fprintf(out, '\n');

% Number of points per cluster.
P_npts_per_clust = ...
SC.plotCombined(arrayfun(@(i) rdata{i}.collected.n_pts, 1 : n_files, ...
                         'UniformOutput', false),                    ...
                1, 'localizations per cluster', econd, 'npts_per_clust');
fprintf(out, 'P_n points per cluster =\n');
fprintf(out, fmt, P_npts_per_clust);
fprintf(out, '\n');

% Nearest neighbour center-to-center cluster distances.
P_nn_dists = ...
SC.plotCombined(arrayfun(@(i) rdata{i}.collected.nn_dists, 1 : n_files, ...
                         'UniformOutput', false),                       ...
                [], 'nn c2c cluster distances (nm)', econd, 'nn_dists');
fprintf(out, 'P_nn c2c cluster distances =\n');
fprintf(out, fmt, P_nn_dists);
fprintf(out, '\n');

% Mean nearest neighbor center-to-center cluster distances per ROI.
P_nn_dists_ROI = ...
SC.plotCombined(arrayfun(@(i) rdata{i}.collected.nn_dists_ROI, 1 : n_files, ...
                         'UniformOutput', false),                           ...
                1, 'mean nn c2c cluster distances per ROI (nm)',            ...
                econd, 'nn_dists_ROI');
fprintf(out, 'P_nn mean c2c cluster distances per ROI =\n');
fprintf(out, fmt, P_nn_dists_ROI);
fprintf(out, '\n');

% Cluster sizes (equivalent radii corresponding to the cluster areas).
P_equiv_radii = ...
SC.plotCombined(arrayfun(@(i) rdata{i}.collected.equiv_radii, 1 : n_files, ...
                         'UniformOutput', false),                          ...
                1, 'cluster equiv radii (nm)', econd, 'equiv_radii');
fprintf(out, 'P_cluster equiv radii =\n');
fprintf(out, fmt, P_equiv_radii);
fprintf(out, '\n');

% Nearest neighbor localization distances.
if CSVglobal
   SC.CSV = true;
end
P_nn_loc_dists = ...
SC.plotCombined(arrayfun(@(i) rdata{i}.collected.nn_loc_dists, 1 : n_files, ...
                         'UniformOutput', false),                           ...
                1, 'nn localization distances (nm)', econd, 'nn_local');
SC.CSV = false;
fprintf(out, 'P_nn localization distances =\n');
fprintf(out, fmt, P_nn_loc_dists);
fprintf(out, '\n');

% Nearest neighbor localization distances within clusters.
P_nn_withinC = ...
SC.plotCombined(arrayfun(@(i) rdata{i}.collected.nn_within_clust, 1:n_files, ...
                         'UniformOutput', false),                            ...
                1, 'nn localization distances within clusters (nm)',         ...
                econd, 'nn_withinC');
fprintf(out, 'P_nn localization distances within clusters =\n');
fprintf(out, fmt, P_nn_withinC);
fprintf(out, '\n');

% ------------------------------------------------
% How to add additional plots using rdata.results:
% ------------------------------------------------
% Clustered fraction per ROI.
clustered_frac = cell(1, n_files);
for j = 1 : n_files
   clustered_frac{j} = [];
   for i = 1 : numel(rdata{j}.results)
      clustered_frac{j} =    ...
         [clustered_frac{j}, ...
          rdata{j}.results{i}.n_clustered / rdata{j}.results{i}.n_points];
   end
end
P_clustered_frac = ...
SC.plotCombined(clustered_frac,                                      ...
                0.05, 'fraction of clustered localizations per ROI', ...
                econd, 'clustered_frac_ROI');
fprintf(out, 'P_clustered fraction per ROI =\n');
fprintf(out, fmt, P_clustered_frac);
fprintf(out, '\n');

% Cluster areas.
areas = cell(1, n_files);
for j = 1 : n_files
   areas{j} = [];
   for i = 1 : numel(rdata{j}.results)
      areas{j} = [areas{j}, rdata{j}.results{i}.areas];
   end
end
P_areas = ...
SC.plotCombined(areas, ...
                5, 'cluster areas (nm^2)', ...
                econd, 'areas');
fprintf(out, 'P_cluster areas =\n');
fprintf(out, fmt, P_areas);
fprintf(out, '\n');

% Cluster areas per ROI.
areas_ROI = cell(1, n_files);
for j = 1 : n_files
   areas_ROI{j} = [];
   for i = 1 : numel(rdata{j}.results)
      areas_ROI{j} = [areas_ROI{j}, mean(rdata{j}.results{i}.areas)];
   end
end
P_areas_ROI = ...
SC.plotCombined(areas_ROI, ...
                5, 'cluster areas per ROI (nm^2)', ...
                econd, 'areas_ROI');
fprintf(out, 'P_cluster areas per ROI =\n');
fprintf(out, fmt, P_areas_ROI);
fprintf(out, '\n');

% Cluster compactness.
compactness = cell(1, n_files);
for j = 1 : n_files
   compactness{j} = [];
   for i = 1 : numel(rdata{j}.results)
      compactness{j} = [compactness{j}, rdata{j}.results{i}.compactness];
   end
end
P_compactness = ...
SC.plotCombined(compactness,                        ...
                0.05, 'cluster compactness', econd, 'compactness');
fprintf(out, 'P_cluster compactness =\n');
fprintf(out, fmt, P_compactness);
fprintf(out, '\n');

% Cluster compactness per ROI.
compactness_ROI = cell(1, n_files);
for j = 1 : n_files
   compactness_ROI{j} = [];
   for i = 1 : numel(rdata{j}.results)
      compactness_ROI{j} = [compactness_ROI{j}, ...
                            mean(rdata{j}.results{i}.compactness)];
   end
end
P_compactness_ROI = ...
SC.plotCombined(compactness_ROI, ...
                0.05, 'cluster compactness per ROI', econd, 'compactness_ROI');
fprintf(out, 'P_cluster compactness per ROI =\n');
fprintf(out, fmt, P_compactness_ROI);
fprintf(out, '\n');

% Number of points per cluster per ROI.
n_pts_ROI = cell(1, n_files);
for j = 1 : n_files
   n_pts_ROI{j} = [];
   for i = 1 : numel(rdata{j}.results)
      n_pts_ROI{j} = [n_pts_ROI{j}, mean(rdata{j}.results{i}.n_pts)];
   end
end
P_npts_per_clust_ROI = ...
SC.plotCombined(n_pts_ROI, 1, 'localizations per cluster per ROI', ...
                econd, 'npts_per_clust_ROI');
fprintf(out, 'P_n points per cluster per ROI =\n');
fprintf(out, fmt, P_npts_per_clust_ROI);
fprintf(out, '\n');

% Number of points per cluster per cluster area.
npts_per_clust_per_area = cell(1, n_files);
for j = 1 : n_files
   n_pts_per_clust_per_area{j} = [];
   for i = 1 : numel(rdata{j}.results)
      n_pts_per_clust_per_area{j} = [n_pts_per_clust_per_area{j}, ...
         rdata{j}.results{i}.n_pts ./ rdata{j}.results{i}.areas];
   end
end
P_n_pts_per_clust_per_area = ...
SC.plotCombined(n_pts_per_clust_per_area, 1,             ...
                'cluster points / area (1/nm^2)', econd, ...
                'npts_per_clust_per_area');
fprintf(out, 'P_cluster points / area =\n');
fprintf(out, fmt, P_n_pts_per_clust_per_area);
fprintf(out, '\n');

% Number of points per cluster per cluster area per ROI.
npts_per_clust_per_area_ROI = cell(1, n_files);
for j = 1 : n_files
   n_pts_per_clust_per_area_ROI{j} = [];
   for i = 1 : numel(rdata{j}.results)
      n_pts_per_clust_per_area_ROI{j} = [n_pts_per_clust_per_area_ROI{j}, ...
         mean(rdata{j}.results{i}.n_pts ./ rdata{j}.results{i}.areas)];
   end
end
P_n_pts_per_clust_per_area_ROI = ...
SC.plotCombined(n_pts_per_clust_per_area_ROI, 1,                 ...
                'cluster points / area (1/nm^2) per ROI', econd, ...
                'npts_per_clust_per_area_ROI');
fprintf(out, 'P_cluster points / area per ROI =\n');
fprintf(out, fmt, P_n_pts_per_clust_per_area_ROI);
fprintf(out, '\n');

% Cluster sizes (equivalent radii corresponding to the cluster areas) per ROI.
equiv_radii_ROI = cell(1, n_files);
for j = 1 : n_files
   equiv_radii_ROI{j} = [];
   for i = 1 : numel(rdata{j}.results)
      equiv_radii_ROI{j} = [equiv_radii_ROI{j}, ...
                            mean(rdata{j}.results{i}.equiv_radii)];
   end
end
%SC.CSV = true;
P_equiv_radii_ROI = ...
SC.plotCombined(equiv_radii_ROI, 1,                        ...
                'cluster equiv radii per ROI (nm)', econd, ...
                'equiv_radii_ROI');
%SC.CSV = false;
fprintf(out, 'P_cluster equiv radii per ROI =\n');
fprintf(out, fmt, P_equiv_radii_ROI);
fprintf(out, '\n');

% Nearest neighbor localization distances per ROI.
nn_loc_dists_ROI = cell(1, n_files);
for j = 1 : n_files
   nn_loc_dists_ROI{j} = [];
   for i = 1 : numel(rdata{j}.results)
      XY_ROI = rdata{j}.results{i}.XY;
      nn_loc_dists_ROI{j} = [nn_loc_dists_ROI{j}, ...
                             mean(c.nn_distances(XY_ROI))];
   end
end
if CSVglobal
   SC.CSV = true;
end
P_nn_loc_dists_ROI = ...
SC.plotCombined(nn_loc_dists_ROI, 1,                             ...
                'nn localization distances per ROI (nm)', econd, ...
                'nn_local_ROI');
SC.CSV = false;
fprintf(out, 'P_nn localization distances per ROI =\n');
fprintf(out, fmt, P_nn_loc_dists_ROI);
fprintf(out, '\n');

% Nearest neighbor localization distances within clusters per ROI.
nn_within_clust_ROI = cell(1, n_files);
for j = 1 : n_files
   nn_within_clust_ROI{j} = [];
   for i = 1 : numel(rdata{j}.results)
      nn_within_clust_ROI{j} = [nn_within_clust_ROI{j}, ...
                                mean(rdata{j}.results{i}.nn_within_clust)];
   end
end
%SC.CSV = true;
P_nn_withinC_ROI = ...
SC.plotCombined(nn_within_clust_ROI, 1,                                   ...
                'nn localization distances within clusters per ROI (nm)', ...
                econd, 'nn_withinC_ROI');
%SC.CSV = false;
fprintf(out, 'P_nn localization distances within clusters per ROI =\n');
fprintf(out, fmt, P_nn_withinC_ROI);
fprintf(out, '\n');

% Number of clusters vs number of localizations per ROI.
clust_vs_loc_ROI = cell(1, n_files);
for j = 1 : n_files
   clust_vs_loc_ROI{j} = [];
   for i = 1 : numel(rdata{j}.results)
      clust_vs_loc_ROI{j} = [clust_vs_loc_ROI{j}, ...
         rdata{j}.results{i}.nC / rdata{j}.results{i}.n_points];
   end
end
%SC.CSV = true;
P_clust_vs_loc_ROI = ...
SC.plotCombined(clust_vs_loc_ROI, 1,                    ...
                'cluster / localization ratio per ROI', ...
                econd, 'clust_vs_loc_ROI');
%SC.CSV = false;
fprintf(out, 'P_cluster / localization ratio per ROI =\n');
fprintf(out, fmt, P_clust_vs_loc_ROI);
fprintf(out, '\n');

fclose(out);

fprintf('Done %s.\n', base_name);

end
