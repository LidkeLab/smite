function [nC, C] = cluster_voronoi(i_rho, self_nbrs, epsilon, minPts, XY)
% Taking the density indices i_rho that identify points to be clustered,
% generate the clusters C (their number given by nC).
%
% INPUTS:
%    i_rho       indices of all the cells exceeding the density criterion
%    self_nbrs   for each cell, the indices for itself and all its neighbors
%    epsilon     epsilon or cutoff distance (nm)
%    minPts      minimum number of points allowed in a cluster
%    XY          point coordinates [N x 2] (nm)
%
% OUTPUTs:
%    nC          number of clusters found
%    C           cell array of XY indices forming each cluster [nC x 1] (nm)

% Created by
%    Michael J. Wester (2019)

   nC = 0;
   C = {};
   while ~isempty(i_rho)
      % There is at least one more point left to be clustered, so create a new
      % cluster and stuff the point in it while deleting it off the list of
      % points remaining to be clustered.
      C_nC = i_rho(1); 
      i_rho(1) = [];
      % Find the point's neighbors and see if any of them are on the list of
      % points remaining to be clustered.
      lo = 1;
      hi = 1;
      nbrs = self_nbrs{C_nC(lo:hi)};
      l = intersect(nbrs, i_rho);
      while ~isempty(l)
         % If there is an overlap between the neighbors of the points just
         % added to the current cluster and those remaining to be clustered,
         % stuff the overlap in the current cluster, delete them off the list
         % of points remaining to be clustered, and then compute the overlap
         % between the neighbors of these new cluster points and those
         % remaining.
         C_nC = [C_nC, l];
         i_rho = setdiff(i_rho, l);
         lo = hi + 1;
         hi = hi + length(l);
         % Find the unique neighbors of the new points just added to the
         % cluster (and themselves).
         nbrs = unique([ self_nbrs{C_nC(lo:hi)} ]);
         l = intersect(nbrs, i_rho);
      end
      % Eliminate clusters smaller than minPts.
      if length(C_nC) >= minPts
         nC = nC + 1;  
         C{nC} = sort(C_nC);
      end
   end

   % If epsilon > 0, apply further restrictions on the clusters found above.
   % For each Voronoi cluster, apply a secondary clustering algorithm that
   % separates points based on epsilon.  This will, in general, separate the
   % Voronoi clusters into smaller (or same size) clusters and isolated points.
   % Collect together the newly separated clusters and return these.
   if epsilon > 0
      c = smi_cluster.Clustering();
      Algorithm = 'Hierarchal';
      nB = 0;
      B  = {};
      for i = 1 : nC
         xy = XY(C{i}, :);
         [nCC, CC, ~, ~] = c.cluster(Algorithm, xy, epsilon, minPts);
         for j = 1 : nCC
            nB = nB + 1;
            % CC{j} are the point indices of a new cluster, so map these back
            % into the point indices of the original cluster C{i}.
            B{nB} = C{i}(CC{j});
         end
      end

      nC = nB;
      C  = B;
   end

end
