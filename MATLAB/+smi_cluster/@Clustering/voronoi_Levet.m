function [area, rho, nC, C] = ...
   voronoi_Levet(obj, xy, alpha, epsilon, minPts, algorithm)
%voronoi_Levet implements Voronoi diagram based clustering.  See citation.
%
% INPUTS:
%    obj         various properties used by the algorithms
%       Plotting     produce Voronoi plots
%       PtIDs        label the points in the plots
%       ShrinkFactor 'boundary' compactness factor in the range 0 - 1
%    xy          point coordinates [N x 2] (nm)
%    alpha       ratio of local density / overall density for a point's Voronoi
%                region to be considered sufficiently dense for clustering
%                purposes
%    epsilon     distance constraint for clustering; no constraint if <= 0 (nm)
%    minPts      minimum number of points needed to form a cluster
%    algorithm   cell array of algorithms to apply:
%                1   [0] calculations consider Voronoi regions only
%                2   [1] calculations consider Voronoi regions and their
%                    adjacent neighbors [DEFAULT]
%                3   [1M] consider the median density of each cell and its
%                    neighbors
%
% OUTPUTS:
%    area   cell array of Voronoi cell areas per point
%    rho    cell array of Voronoi cell densities per point
%    nC     cell array of number of clusters discovered
%    C      cell array of indices of points in each cluster
%
% CITATION:
%    Florian Levet, Eric Hosy, Adel Kechkar, Corey Butler, Anne Beghin, Daniel
%    Choquet and Jean-Baptiste Sibarita, ``SR-Tesseler: a method to segment and
%    quantify localization-based super-resolution microscopy data'', _Nature
%    Methods_, 2015 (DOI:10.1038/NMETH.3579).

% Created by
%    Michael J. Wester (2019)

   if ~exist('algorithm', 'var')
      algorithm = 2;
   end

   if any(algorithm < 1 | algorithm > 3)
      error('Invalid algorithm selected: %d', algorithm);
   end

   ShrinkFactor = obj.ShrinkFactor;

   XY = double(xy);
   dim = size(XY, 2);

   X = XY(:, 1);
   Y = XY(:, 2);
   n_XY = length(X);
   if dim == 2
      area_XY = (max(X) - min(X)) * (max(Y) - min(Y));
   else
      Z = XY(:, 3);
      area_XY = (max(X) - min(X)) * (max(Y) - min(Y)) * (max(Z) - min(Z));
   end
   rho_XY = n_XY / area_XY;

   [v, c] = voronoin(XY);
   n_v = size(v, 1);
   n_c = length(c);

   % Find what Voronoi cells (c) contain each vertex (v).
   v2c = cell(1, n_v);  
   for i = 1 : n_c
      c_i = c{i};
      for j = 1 : length(c_i)
         v2c{c_i(j)} = [v2c{c_i(j)}, i];
      end
   end

   % For each cell, collect together the indices for itself and all its
   % neighbors.
   self_nbrs = cell(1, n_c);
   for i = 1 : n_c
      c_i = c{i};
      for j = 1 : length(c_i)
         k = v2c{c_i(j)};
         self_nbrs{i} = [self_nbrs{i}, k];
      end
      self_nbrs{i} = unique(self_nbrs{i});
   end

   % Level 0: consider each cell individually.
   area_0 = zeros(1, n_XY);
   rho_0  = zeros(1, n_XY);
   for i = 1 : n_c
      c_i = c{i};
      if all(c_i ~= 1)
         if dim == 2
            area_0(i) = polyarea(v(c_i, 1), v(c_i, 2));
         else
            [~, area_0(i)] = boundary(v(c_i, 1), v(c_i, 2), v(c_i, 3), ...
                                      ShrinkFactor);
         end
         rho_0(i) = 1 / area_0(i);
      else
         area_0(i) = Inf;
         rho_0(i) = 0;
      end
   end

   if any(algorithm == 1)
      [~, i_rho_0] = find(rho_0 >= alpha*rho_XY);
      [nC_0, C_0] = smi_cluster.Clustering.cluster_voronoi( ...
                       i_rho_0, self_nbrs, epsilon, minPts, XY);

      if obj.Plotting
         if dim == 2
            smi_cluster.Clustering.plot_voronoi(X, Y, v, c, ...
               rho_0 ./ rho_XY, 'rho_0 / rho_a', i_rho_0, obj.PtIDs);
         else
            smi_cluster.Clustering.plot_voronoi3(X, Y, Z, v, c, ...
               rho_0 ./ rho_XY, 'rho_0 / rho_a', i_rho_0, obj.PtIDs);
         end
      end

      area{1} = area_0;
      rho{1}  = rho_0;
      nC{1}   = nC_0;
      C{1}    = C_0;
   end

   if any(algorithm == 2)
      % Level 1: consider each cell and its neighbors.
      area_1 = zeros(1, n_XY);
      rho_1  = zeros(1, n_XY);
      for i = 1 : n_c
         c_i = c{i};
         if all(c_i ~= 1)
            % Cells with vertex 1, which is the point at Infinity, are
            % excluded.
            n_self_nbrs = length(self_nbrs{i});
            for j = 1 : n_self_nbrs
               area_1(i) = area_1(i) + area_0(self_nbrs{i}(j));
            end
            rho_1(i) = n_self_nbrs / area_1(i);
         else
            area_1(i) = Inf;
            rho_1(i) = 0;
         end
      end

      [~, i_rho_1] = find(rho_1 >= alpha*rho_XY);
      [nC_1, C_1] = smi_cluster.Clustering.cluster_voronoi( ...
                       i_rho_1, self_nbrs, epsilon, minPts, XY);

      if obj.Plotting
         if dim == 2
            smi_cluster.Clustering.plot_voronoi(X, Y, v, c, ...
               rho_1 ./ rho_XY, 'rho_1 / rho_a', i_rho_1, obj.PtIDs);
         else
            smi_cluster.Clustering.plot_voronoi3(X, Y, Z, v, c, ...
               rho_1 ./ rho_XY, 'rho_1 / rho_a', i_rho_1, obj.PtIDs);
         end
      end

      area{2} = area_1;
      rho{2}  = rho_1;
      nC{2}   = nC_1;
      C{2}    = C_1;

      %area_1A = zeros(1, n_XY);
      %rho_1A  = zeros(1, n_XY);
      %for i = 1 : n_c
      %   c_i = c{i};
      %   if all(c_i ~= 1)
      %      n_self_nbrs = 0;
      %      for j = 1 : n_c
      %         if ~isempty(intersect(c{i}, c{j}))
      %            n_self_nbrs = n_self_nbrs + 1;
      %            area_1A(i) = area_1A(i) + area_0(j);
      %         end
      %      end
      %      rho_1A(i) = n_self_nbrs / area_1A(i);
      %   else
      %      area_1A(i) = Inf;
      %      rho_1A(i) = 0;
      %   end
      %end
      %plot_voronoi(X, Y, v, c, rho_1A);
   end

   if any(algorithm == 3)
      % Level 1M: consider the median density of each cell and its neighbors.
      area_1M = area_0;
      rho_1M  = zeros(1, n_XY);
      for i = 1 : n_c
         rho_1M(i) = median(rho_0(self_nbrs{i}));
      end

      [~, i_rho_1M] = find(rho_1M >= alpha*rho_XY);
      [nC_1M, C_1M] = smi_cluster.Clustering.cluster_voronoi( ...
                         i_rho_1M, self_nbrs, epsilon, minPts, XY);

      if obj.Plotting
         if dim == 2
            smi_cluster.Clustering.plot_voronoi(X, Y, v, c, ...
               rho_1M ./ rho_XY, 'rho_{1M} / rho_a', i_rho_1M, obj.PtIDs);
         else
            smi_cluster.Clustering.plot_voronoi3(X, Y, Z, v, c, ...
               rho_1M ./ rho_XY, 'rho_{1M} / rho_a', i_rho_1M, obj.PtIDs);
         end
      end

      area{3} = area_1M;
      rho{3}  = rho_1M;
      nC{3}   = nC_1M;
      C{3}    = C_1M;
   end

end
