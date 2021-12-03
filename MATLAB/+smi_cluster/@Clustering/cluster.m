function [nC, C, centers, ptsI] = ...
   cluster(obj, algorithm, SMD, E, minPts)
% Main interface to the clustering algorithms described below:
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
%
% INPUTS:
%    obj         various properties used by the algorithms
%                --- Properties used by voronoi_Levet: ---
%       Alpha       [2] ratio of local density / overall density for a
%                   point's Voronoi region to be considered sufficiently
%                   dense for clustering purposes
%       Valgorithm  [2] Voronoi algorithm to apply:
%                      1   [0] calculations consider Voronoi regions only
%                      2   [1] calculations consider Voronoi regions and
%                              their adjacent neighbors
%                      3   [1M] consider the median density of each cell
%                               and its neighbors
%       Plotting    [false] produce Voronoi plots
%       PtIDs       [false] label the points in the plots
%                --- Properties used by H-SET: ---
%       PixelSize   [100]      pixel size (nm)
%       Sigma_Reg   [[10, 10]] x, y registration error (nm) 
%    algorithm   one of the clustering algorithms below:
%                DBSCAN_Daszykowski or DBSCAN
%                DBSCAN_Daszykowski_noE (computes its own value for E)
%                Hierarchical or Hierarchal
%                           (Matlab's hierarchical clustering algorithm)
%                Voronoi    (Florian Levet et al's Voronoi based algorithm)
%                H-SET      (clustering implied by the H-SET algorithm)
%    SMD         (x, y) coordinates in the format
%                    (1) SMD structure: SMD.X, SMD.Y with optional SMD.X_SE,
%                        SMD.Y_SE (only needed for H-SET clustering) (pixel)
%                    (2) N x 2 or N x 3 array of coordinates (nm)
%    E           epsilon or cutoff distance (nm).  This is the maximum distance
%                between points in a cluster or the minimum distance between
%                points in different clusters
%    minPts      [OPTIONAL] minimum number of points allowed in a cluster
%                (default = 3)
%
% OUTPUTS:
%    nC          number of clusters found
%    C           cell array of XY indices forming each cluster [nC x 1]
%    centers     coordinates of the center of each cluster [nC x n_dim]
%    ptsI        indices of points not found in any cluster

   if ~exist('minPts', 'var')
      minPts = 3;
   end

   % Check for SMD structure; in this situation, assume the units are pixels.
   SMDstruct = false;
   if isfield(SMD, 'X') & isfield(SMD, 'Y')
      SMDstruct = true;
      if isfield(SMD, 'Z') && ~isempty(SMD.Z)
         XY = [SMD.X, SMD.Y, SMD.Z] .* obj.PixelSize;
      else
         XY = [SMD.X, SMD.Y] .* obj.PixelSize;
      end
   % Check for N x 2 or N x 3 matrix of coordinates.
   elseif ismatrix(SMD) && (size(SMD, 2) == 2 || size(SMD, 2) == 3)
      XY = SMD;
   else
      error('Unrecognized coordinate input!');
   end

   if strcmp(algorithm, 'H-SET') && SMDstruct && ...
         (~isfield(SMD, 'X_SE') || ~isfield(SMD, 'Y_SE'))
      error('X_SE/Y_SE needed for H-SET!');
   end

   n_dim = size(XY, 2);

   if size(XY, 1) == 0
      warning('No points to cluster!');

      nC = 0;
      C = cell(1, nC);
      centers = zeros(n_dim, nC);
      ptsI = [];

      return
   end

   switch algorithm
   case {'DBSCAN_Daszykowski', 'DBSCAN'}

      [ptsC, Ctype] = obj.dbscan_Daszykowski(XY, minPts - 1, E);
      nC = max([0, ptsC]);
      C = cell(1, nC);
      centers = zeros(n_dim, nC);
      for j = 1 : nC
         C{j} = find(ptsC == j);
         centers(:, j) = mean(XY(C{j}, :), 1);
      end
      ptsI = find(ptsC <= 0);
      %isolated = XY(ptsI, :);

   % Here, E (epsilon) is calculated internally by the DBSCAN algorithm.  Note
   % that this value of E will be specific to the input data in XY.
   case 'DBSCAN_Daszykowski_noE'

      [ptsC, Ctype, E] = obj.dbscan_Daszykowski(XY, minPts - 1, []);
      nC = max([0, ptsC]);
      C = cell(1, nC);
      centers = zeros(n_dim, nC);
      for j = 1 : nC
         C{j} = find(ptsC == j);
         centers(:, j) = mean(XY(C{j}, :), 1);
      end
      ptsI = find(ptsC <= 0);

   case {'Hierarchical', 'Hierarchal'}

      [C, ptsI] = obj.hierarchal(XY, E, minPts);
      nC = length(C);
      centers = zeros(n_dim, nC);
      for i = 1 : nC
         centers(:, i) = mean(XY(C{i}, :), 1);
      end

   case 'Voronoi'

      %[area, rho, nC, C] = ...
      %   voronoi_Levet(XY, alpha, epsilon, minPts, algorithm);
      [~, ~, nC_a, C_a] = ...
         obj.voronoi_Levet(XY, obj.Alpha, E, minPts, obj.Valgorithm);
      nC = nC_a{obj.Valgorithm};
      C  = C_a{obj.Valgorithm};
      centers = zeros(n_dim, nC);
      for i = 1 : nC
         centers(:, i) = mean(XY(C{i}, :), 1);
      end
      ptsI = setdiff(1 : size(XY, 1), horzcat(C{:}));

   case 'H-SET'

      % XY is assumed to be the same as [SMD.X, SMD.Y] .* obj.PixelSize
      [nC, C] = obj.cluster_HSET(SMD, minPts);
      centers = zeros(n_dim, nC);
      for i = 1 : nC
         centers(:, i) = mean(XY(C{i}, :), 1);
      end
      ptsI = setdiff(1 : size(XY, 1), horzcat(C{:}));

   otherwise

      error('Unknown algorithm: %s\n', algorithm);

   end

end
