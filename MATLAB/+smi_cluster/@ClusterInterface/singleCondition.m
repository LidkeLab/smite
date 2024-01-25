function singleCondition(pathname, files, algorithm_range, E_range, ...
                         minPts_range, Pixel2nm, base_name, A_ROI,  ...
                         doHopkins, doSigmaActual, Alpha)
% ---------- Statistics for a single condition
% Perform cluster analysis for comparison of experimental conditions.  This is
% useful for performing parameter studies, especially over a range of E's
% (epsilons) and minPts' (N's).
%
% INPUTS:
%    pathname          path to where the files are located
%    files             _ROIs.mat files to collect ROI info from
%    algorithm_range   clustering algorithms to try (cell array of strings):
%                          {'DBSCAN', 'Hierarchical', 'Voronoi'}
%    E_range           minimum distances between clusters or maximum distance
%                      between points in a cluster (nm)
%    minPts_range      minimum numbers of points in a cluster
%    Pixel2nm          conversion factor from pixels to nm
%    base_name         analysis identification string
%    A_ROI             ROI area (nm^2)
%    doHopkins         [DEFAULT = true] Hopkins' test can be time consuming for
%                      dense ROIs
%    doSigmaActual     [DEFAULT = true] set to false for very dense ROIs to
%                      avoid crashes due to lack of memory
%    Alpha             [DEFAULT = 2] ratio of local density / overall density
%                      for Voronoi clustering
%                     
%
% OUTPUTS:
%    A variety of plots are produced and are placed in a subdirectory of
%    pathname, identified by the algorithm, minPts and E value used for
%    clustering, e.g., DBSCAN_N=3,E=30.  Plots are produced for each ROI of
%    the localizations and clusters detected [boundary colors are random]
%    (*_ROI#?_clusters), and a plot comparing the nearest neighbor histogram
%    of the localizations with a random distribution (_ROI#?_NN_loc_PDF).  A
%    variety of statistics are also computed and plotted (with prefix ALL_)
%    over the condition, and a .mat file is saved with a name like
%       base_nam#algorithm_N=*,E=*_results.mat
%    which is just the naming convention needed for comparing results of
%    multiple conditions by combinedStatistics1/2).

% Created by
%    Michael J. Wester (2022)

   if ~exist('doHopkins', 'var')
      doHopkins = true;
   end
   if ~exist('doSigmaActual', 'var')
      doSigmaActual = true;
   end
   if ~exist('Alpha', 'var')
      Alpha = [];
   end

   c  = smi_cluster.Clustering();
   SC = smi_cluster.StatisticsClustering();
   PC = smi_cluster.PairCorrelation();

   base_name_save = base_name;
   if ~startsWith(base_name, 'ALL')
      base_name = ['ALL_', base_name];
   end
   SC.BaseName = base_name;
   PC.BaseName = base_name;

   % Set to false for very dense ROIs to avoid crashes due to lack of memory.
   c.DoSigmaActual = doSigmaActual;
   if ~isempty(Alpha)
      c.Alpha = Alpha;
   else
      c.Alpha = 2;
   end
   %c.Alpha = 800;   c.Plotting = true;

   LoS = 0.01;   % level of significance for H-SET
   c.LoS = LoS;
   User_Sigma_Reg = [10, 10];   % specified image registration error for H-SET

   % Cluster boundary is halfway between a convex hull and fully concave.
   shrinkFactor = 0.5;
   c.ShrinkFactor = shrinkFactor;

   for E = E_range
   for minPts = minPts_range
   for alg = algorithm_range
   algorithm = alg{1};

   % Make a subdirectory for this particular analysis.
   switch algorithm
   case 'Voronoi'
      Run = sprintf('%s_N=%d,Alpha=%g', algorithm, minPts, c.Alpha);
      Run = regexprep(Run, '\.', '_');
   case 'H-SET'
      %Run = sprintf('%s_N=%d,LoS=%0.2f', algorithm, minPts, c.LoS);
      Run = sprintf('%s_N=%d,LoS=%g,SReg=%d', ...
                    algorithm, minPts, c.LoS, User_Sigma_Reg(1));
   otherwise
      Run = sprintf('%s_N=%d,E=%d', algorithm, minPts, E);
   end
   Aresults_dir = fullfile(pathname, Run);
   % Create Aresults_dir if it does not already exist.
   if ~isfolder(Aresults_dir)
      mkdir(Aresults_dir);
   end

   n_files = numel(files);
   if n_files == 0
      error('No files to analyze in %s!', pathname);
   end

   c.ResultsDir  = Aresults_dir;
   SC.ResultsDir = Aresults_dir;
   PC.ResultsDir = Aresults_dir;

   % Pair/auto-correlation properties.
   PC.Fig_ext = 'png';   % figure extension
   PC.Rmax_axis = 500;   % Sets plotting limit if > 0
   % Histogram bin size for pairwise correlation---this is the number of pixels
   % per bin over which correlation statistics are collected.
   PC.HistBinSize = 5;
   PC.Verbose = 0;       % verbose output and extra saved .mat files

   % Cluster all ROIs in all images, then collect statistics.
   n_ROIs = [];
   RoI = {};
   nC = 0;
   n_clusts = [];
   c_density = [];
   n_points = [];
   n_pts = [];
   compactness = [];
   nn_dists = [];
   nn_dists_ROI = [];
   equiv_radii = [];
   nn_loc_dists = [];
   nn_within_clust = [];
   clear results ROIs

   n = 0;   % count the total number of ROIs
   for j = 1 : n_files
      data = load(fullfile(pathname, files{j}));
      short = regexprep(files{j}, '_ROIs.mat$', '');
      n_ROIs = [n_ROIs, data.n_ROIs];
      RoI{j} = data.RoI;
      for i = 1 : data.n_ROIs
         nPts = numel(data.RoI{i}.X{1});
         % Skip completely empty ROIs!
         if nPts == 0
            continue;
         end
         n = n + 1;
         XY_ROI = [data.RoI{i}.X{1}, data.RoI{i}.Y{1}];
         %ROIs{n}.XY1 = XY_ROI;
         ROIs{n}.X = { data.RoI{i}.X{1} };
         ROIs{n}.Y = { data.RoI{i}.Y{1} };
         ROIs{n}.ROI = data.RoI{i}.ROI;
         fprintf('%s file %d ROI %d number of points = %d\n', Run, j, i, nPts);

         if strcmp(algorithm, 'H-SET')
            SMR = {};
            SMR.X = data.RoI{i}.X{1} ./ Pixel2nm;
            SMR.Y = data.RoI{i}.Y{1} ./ Pixel2nm;
            SMR.X_SE = data.RoI{i}.X_STD{1} ./ Pixel2nm;
            SMR.Y_SE = data.RoI{i}.Y_STD{1} ./ Pixel2nm;
            if exist('User_Sigma_Reg', 'var') && ~isempty('User_Sigma_Reg')
               SMR.Sigma_Reg = User_Sigma_Reg;
            else
               SMR.Sigma_Reg = data.Sigma_Reg;
            end
            [nC, C, centers, ptsI] = ...
               c.cluster(algorithm, XY_ROI, E, minPts, SMR);
         else
            [nC, C, centers, ptsI] = ...
               c.cluster(algorithm, XY_ROI, E, minPts);
         end
         fprintf('%s file %d ROI %d number of clusters = %d\n', Run, j, i, nC);
         results{n} = c.clusterStats(XY_ROI, C, centers);

         nC = nC + results{n}.nC;
         n_clusts     = [n_clusts,     results{n}.nC];
         c_density    = [c_density,    results{n}.nC / A_ROI];
         n_points     = [n_points,     results{n}.n_points];
         n_pts        = [n_pts,        results{n}.n_pts];
         compactness  = [compactness,  results{n}.compactness];
         nn_dists     = [nn_dists,     results{n}.min_c2c_dists];
         nn_dists_ROI = [nn_dists_ROI, mean(results{n}.min_c2c_dists)];
         equiv_radii  = [equiv_radii,  results{n}.equiv_radii];
         nn_loc_dists = [nn_loc_dists, c.nn_distances(XY_ROI)];
         nn_within_clust = [nn_within_clust, results{n}.nn_within_clust];

         clusterFig = c.plotClusters(XY_ROI, C, centers, ptsI, algorithm);
%        showm(clusterFig);   % comment out if not displaying to the screen
         filename = fullfile(Aresults_dir, ...
                             sprintf('%s_ROI#%d_clusters', short, i));
         saveas(clusterFig, filename, 'fig');
         %saveas(clusterFig, filename, 'png');
         try
            print(clusterFig, '-r300', filename, '-dpng');
         catch ME
            fprintf('### PROBLEM with %s ###\n', filename);
            fprintf('%s\n', ME.identifier);
            fprintf('%s\n', ME.message); 
         end

         if strcmp(algorithm, 'H-SET')
            XY_SE_ROI = [data.RoI{i}.X_STD{1}, data.RoI{i}.Y_STD{1}];
            clusterFig = c.plotClustersSE(XY_ROI, XY_SE_ROI, C, centers, ...
                                          ptsI, algorithm);
%           showm(clusterFig);   % comment out if not displaying to the screen
            filename = fullfile(Aresults_dir, ...
                                sprintf('%s_ROI#%d_SEclusters', short, i));
            %saveas(clusterFig, filename, 'fig');
            %saveas(clusterFig, filename, 'png');
            if nPts < 25000
               try
                  print(clusterFig, '-r600', filename, '-dpng');
               catch ME
                  fprintf('### clusterFig PROBLEM with %s ###\n', filename);
                  fprintf('%s\n', ME.identifier);
                  fprintf('%s\n', ME.message); 
               end
            end
         end

         c.nn_ROIrandom(XY_ROI, A_ROI, ...
                        sprintf('%s_ROI#%d localizations', short, i));
         filename = fullfile(Aresults_dir, ...
                             sprintf('%s_ROI#%d_NN_loc_PDF', short, i));
         try
            print(gcf, '-r300', filename, '-dpng');
         catch ME
            fprintf('### nn_ROIrandom PROBLEM with %s ###\n', filename);
            fprintf('%s\n', ME.identifier);
            fprintf('%s\n', ME.message); 
         end

         close
      end % i
   end % j

   collected.A_ROI = A_ROI;
   collected.nC = nC;
   collected.n_clusts = n_clusts;
   collected.c_density = c_density;
   collected.n_points = n_points;
   collected.n_pts = n_pts;
   collected.nn_dists = nn_dists;
   collected.nn_dists_ROI = nn_dists_ROI;
   collected.equiv_radii = equiv_radii;
   collected.nn_loc_dists = nn_loc_dists;
   collected.nn_within_clust = nn_within_clust;

   % Hopkins' statistic per ROI.
   if doHopkins
      H_ROI = SC.hopkins_ROIcombined(n, ROIs);
      collected.H_ROI = H_ROI;
   end

   % Save statistical results in case they are needed later.
%  save(fullfile(Aresults_dir, [base_name, '_results.mat']), ...
%       'files', 'n_ROIs', 'RoI', 'Pixel2nm', 'results', 'collected');
   save(fullfile(Aresults_dir, [base_name_save, '#', Run, '_results.mat']), ...
        'files', 'n_ROIs', 'RoI', 'Pixel2nm', 'results', 'collected');

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
   SC.PlotDo = 'pCSx';
   % Red mean, green median (2 only mean, 3 only median) for PlotSpread plots.
   %SC.ShowMM = 1;
   % Options for CDF2 plots are: 'plot', 'semilogx', 'semilogy', 'loglog'.
   %SC.LinLog = 'semilogx';

   econd = {};

   % Number of clusters per ROI.
   SC.plotCombined({n_clusts}, 1, '# of clusters per ROI', econd, ...
                  '_n_clusts_ROI');

   % Density of clusters per ROI.
   SC.plotCombined({c_density}, 1/(A_ROI/1000^2),                     ...
                   'density of clusters per ROI (\mu m^{-2})', econd, ...
                   '_c_density_ROI');

   % Number of localizations per ROI.
   SC.plotCombined({n_points}, [], '# of localizations per ROI', econd, ...
                   '_n_local_ROI');

   % Number of points per cluster.
   SC.plotCombined({n_pts}, 1, 'localizations per cluster', econd, ...
                   '_npts_per_clust');

   % Cluster compactness.
   SC.plotCombined({compactness}, 1, 'cluster compactness', econd, ...
                   '_compactness');

   % Nearest neighbour center-to-center cluster distances.
   SC.plotCombined({nn_dists}, [], 'nn c2c cluster distances (nm)', econd, ...
                   '_nn_dists');

   % Cluster sizes (equivalent radii corresponding to the cluster areas).
   SC.plotCombined({equiv_radii}, 1, 'cluster equiv radii (nm)', econd, ...
                   '_equiv_radii');

   % Nearest neighbor localization distances.
   SC.plotCombined({nn_loc_dists}, 1, 'nn localization distances (nm)', ...
                   econd, '_nn_local');

   % Nearest neighbor localization distances within clusters.
   SC.plotCombined({nn_within_clust}, 1,                             ...
                   'nn localization distances within clusters (nm)', ...
                   econd, '_nn_withinC');

   % Auto-correlation
   fprintf('ROI combined\n\n');
   txt = sprintf('%s_RC', base_name);
   resultsRC_pcc = PC.pair_correlation_ROIcombined(1, n, ROIs, 1);

   end % alg
   end % minPts
   end % E

   clear User_Sigma_Reg

   fprintf('Done %s_%s.\n', base_name, Run);

end
