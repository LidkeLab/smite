### +smi_cluster/@StatisticsClustering

StatisticsClustering contains statistics and plots for detecting clustering.
Data can be 2D or 3D.

See [MATLAB/examples/Example_StatisticsClustering.m](../../examples/Example_StatisticsClustering.m)
for 2D and 3D examples.

The main routines here are:
```
   hopkins              hopkins_ROIcombined           Hopkins' statistic
   ripley               ripley_ROIcombined            Ripley's statistics
   bivariateRipley      bivariateRipley_ROIcombined   Bivariate Ripley's
   pairwiseDist                                       Pairwise distances
   pairwiseMutualDist                                 Piarwise mutual dist.
   plotCombined
      Frequency, CDF, PDF, plotSpread, box and bar plots of an array.
```
The _ROIcombined plots are averaged results over a series of ROIs and
expect n_ROIs and a RoI structure (see smi_helpers.ROITools).

Note that pairwiseMutualDist and bivariateRipley are expecting either a
single argument RoI structure or two SMD structures or two Nx2 (or Nx3)
arrays, while the other statistical functions are expecting an RoI
structure or one SMD structure or one array.

**plotCombined** is a handy way to plot the same data in several different ways
as indicated by the user.

---

properties:
```
   % If ROI is provided, it will be used, otherwise the xy_size will be
   % calculated using the (x_min, x_max, y_min, y_max) computed from the data
   % as
   %    xy_size = min(x_max - x_min, y_max - y_min)
   ROI = [];   % [x_min, x_max, x_max, y_max]   % nm

   Font_props = {'FontSize', 15, 'FontWeight', 'bold'};
   Line_props = {'LineWidth', 3};
   Fig_ext = 'png';
   BaseName = '';         % Descriptive name for the result files
   ResultsDir = '.';      % Directory to store results
   Xlim = [];             % x-axis limits if defined
   Ylim = [];             % y-axis limits if defined

   % Properties used by plotCombined.
   PlotDo = 'fnpcCsSxb';  % Plots to do for plotCombined (fnpcCsSxb)
   LinLog = 'plot';       % Plot type for plotCombined CDF2 plots:
                          %   'plot', 'semilogx', 'semilogy', 'loglog'
   LegendTitle = '';      % Optional legend title
   ShowMM = 1;            % Red mean and green median (2 only mean, 3
                          %    only median) for PlotSpread plots
   CSV = false;           % If true, produce a CSV file of the data

   % Properties below are used by various routines.
   Rate          =  20;   % Sampling rate for statistical functions
   Dendro_cutoff =  50;   % Cluster cutoff for Dendrogram analysis (nm)
   Ripley_cutoff = 200;   % Ripley distance cutoff (nm)
   Confidence    = 2.576; % Bivariate Ripley confidence interval
      % 95% confidence interval: Confidence = 1.96
      % 99% confidence interval: Confidence = 2.576
   Nsims         = 20;    % Simulations to run for bivariateRipley,
                          % pairwiseDist, pairwiseMutualDist

   PixelSize = 100;       % Camera pixel size (nm)

   Verbose = 1;           % verbosity level
```

---

methods:
- **[bivariateRipley](bivariateRipley.m)**:
  Bivariate Ripley's statistic to test the clustering of the points in
  SMD1 and SMD2
- **[bivariateRipley_ROIcombined](bivariateRipley_ROIcombined.m)**:
  Use bivariate Ripley statistics to test the clustering of a series of ROIs
  all of the same size
- **[histogram](histogram.m)**:
  simple histogram function
- **[hopkins](hopkins.m)**:
  Use the Hopkins' statistic to test the clustering of the points in SMD
- **[hopkins_ROIcombined](hopkins_ROIcombined.m)**:
  Use Hopkins' statistics to test the clustering of a series of ROIs
- **[hopkinstat](hopkinstat.m)**:
  Compute a Hopkin's statistic for the 2D particles P
- **[hopkinstat3](hopkinstat3.m)**:
  Compute a Hopkin's statistic for the 3D particles P
- **[pairwiseDist](pairwiseDist.m)**:
  plots pairwise distances and CDFs
- **[pairwiseMutualDist](pairwiseMutualDist.m)**:
  Plot pairwise distances and CDFs between two populations of particles
- **[plotCombined](plotCombined.m)**:
  Combined frequency, CDF, PDF, plotSpread, box and bar plots of the arrays y
- **[ripley](ripley.m)**:
  Use Ripley's statistics to test the clustering of the points in SMD
- **[ripley_ROIcombined](ripley_ROIcombined.m)**:
  Use Ripley's statistics to test the clustering of a series of ROIs
- **[unitTest](unitTest.m)**:
  Test calling StatisticsClustering routines
