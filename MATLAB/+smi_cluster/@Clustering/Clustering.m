classdef Clustering < handle

% Clustering class written by Michael Wester, Keith Lidke, Carolyn Pehlke, Flor
%    Espinoza Hidalgo, Stanly Steinberg and others as noted internally
%    (2/22/2018) <wester@math.unm.edu>
% The New Mexico Center for the Spatiotemporal Modeling of Cell Signaling
% University of New Mexico Health Sciences Center
% Albuquerque, New Mexico, USA   87131
% Copyright (c) 2015-2021 by Michael J. Wester and Keith A. Lidke
%
%    DBSCAN_Daszykowski is the recommended DBSCAN algorithm as it is both fast
%    and stable under coordinate reordering.
%
% Martin Ester, Hans-Peter Kriegel and J\"org Sander and Xiaowei Xu, ``A
% Density-Based Algorithm for Discovering Clusters in Large Spatial Databases
% with Noise'', in _Proceedings of 2nd International Conference on Knowledge
% Discovery and Data Mining (KDD-96)_ edited by Evangelos Simoudis, Jiawei Han
% and Usama M. Fayyad, AAAI Press, 1996, 226--231, ISBN:1-57735-004-9,
% DOI:10.1.1.71.1980.
%
% M. Daszykowski, B. Walczak, D. L. Massart, Looking for Natural Patterns in
% Data. Part 1: Density Based Approach, Chemom. Intell. Lab. Syst. 56 (2001)
% 83-92. 
%
%    Hierarchical is Matlab's hierarchal clustering algorithm with some small
%    additions.
%
%    Voronoi is Florian Levet et al's Voronoi based algorithm.
%
% Florian Levet, Eric Hosy, Adel Kechkar, Corey Butler, Anne Beghin, Daniel
% Choquet and Jean-Baptiste Sibarita, ``SR-Tesseler: a method to segment and
% quantify localization-based super-resolution microscopy data'', _Nature
% Methods_, Volume 12, Number 11, 2015, 1065--1071 (DOI:10.1038/NMETH.3579).
%
%    H-SET is the clustering implied by H-SET, in which the nodes that would be
%    combined in normal H-SET are taken to be clusters here.
%
% Jia Lin, Michael J. Wester, Matthew S. Graus, Keith~A. Lidke and Aaron K.
% Neumann, ``Nanoscopic cell wall architecture of an immunogenic ligand in
% _Candida albicans_ during antifungal drug treatment'', _Molecular
% Biology of the Cell_, Volume 27, Number 6, March 15, 2016, 1002--1014
% (DOI: 10.1091/mbc.E15-06-0355, PMID: 26792838, PMCID: PMC4791122).

% =============================================================================
properties
% =============================================================================

   % --- Generic properties.
   Font_props = {'FontSize', 15, 'FontWeight', 'bold'};
   Fig_ext = 'png';
   ResultsDir = '.'; % Directory to store results.
   ShrinkFactor = 0.5; % 'boundary' compactness factor in the range 0 - 1
   Verbose = 0; % verbosity level
   Xlim = [];   % x-axis limits if defined
   Ylim = [];   % y-axis limits if defined

   % Properties used by clusterStats.
   DoSigmaActual = true; % this calc. can run out of memory for very dense ROIs
                         % (see clusterStats), so make it optional

   % Properties used by voronoi_Levet.
   Alpha      = 2;
      % Ratio of local density / overall density for a point's Voronoi
      % region to be considered sufficiently dense for clustering purposes
   Valgorithm = 2;       % Voronoi algorithm to apply:
      % 1   [0]  calculations consider Voronoi regions only
      % 2   [1]  calculations consider Voronoi regions and their adjacent
      %          neighbors
      % 3   [1M] consider the median density of each cell and its neighbors
   Plotting   = false;   % Produce Voronoi plots
   PtIDs      = false;   % Label the points in the plots

   % --- H-SET properties.
   % H-SET collapse method: 'hierarchal_singlelabel' or 'trivial'.
   Method = 'hierarchal_singlelabel';
   LoS = 0.01;           % level of significance
   PixelSize = 100;      % conversion from pixels to nm
   PlotFigures = true;   % plot various cluster related figures
   Sigma_Reg = [10, 10]; % registration error in x, y (nm)
   Timing = true;        % produce timings for clustering

% =============================================================================
end % properties
% =============================================================================

% =============================================================================
methods
% =============================================================================
% Constructor.  SMF is an optional argument.
   function obj = Clustering(SMF)

      if ~exist('SMF', 'var')
         SMF = smi_core.SingleMoleculeFitting();
      end
      obj.ResultsDir = SMF.Data.ResultsDir;
      obj.PixelSize  = SMF.Data.PixelSize;

   end
% =============================================================================
end % methods
% =============================================================================

% =============================================================================
methods(Static)
% =============================================================================

   [nC, C] = cluster_voronoi(i_rho, self_nbrs, epsilon, minPts, XY)
   [class,type,Eps]=dbscan_Daszykowski(x,k,Eps)
   min_d = edge2edge(hull1, hull2)
   [XY_combined, sigma_combined, nodes_combined] = ...
      hierarchalSingleLabel(XY, sigma, Sigma_Reg, LoS)
   [C, ptsI] = hierarchal(XY, E, minPts)
   [lia] = my_ismemberBuiltinTypes(a,b)
   min_d = nn_distances(xy)
   plot_voronoi(X, Y, v, c, rho, str, dense, ptIDs)
   plot_voronoi3(X, Y, Z, v, c, rho, str, dense, ptIDs)
   [Pvalue, X_Point, Sigma_Point] = singleLabelTest(X, Sigma, Sigma_Reg)
   success = unitTest()
   max_d = vertex2vertex(hull)

% =============================================================================
end % methods(Static)
% =============================================================================
end % classdef
