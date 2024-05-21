function results = clusterStats(obj, SMD, C, centers)
%clusterStats produces statistics on computed clusters.
%
% INPUTS:
%    SMD            (x, y) coordinates of all the points processed (nm) in the
%                   format
%                      (1) SMD structure: SMD.X, SMD.Y with optional SMD.X_SE,
%                          SMD.Y_SE (needed for H-SET clustering) (pixel)
%                      (2) N x 2 or N x 3 array of coordinates (nm)
%    C              cell array of the indices of the points in each cluster
%    centers        array of the coordinates of the center of each cluster
%
% OUTPUTS (contained in results):
%    nC                number of clusters
%    C                 cell array of the indices of the points in each cluster
%    SMD               (x, y) coordinates in the format
%                       (1) SMD structure: SMD.X, SMD.Y with optional SMD.X_SE,
%                           SMD.Y_SE (only needed for H-SET clustering) (pixel)
%                       (2) N x 2 or N x 3 array of coordinates (nm)
%    n_points          total number of points
%    n_clustered       number of points in clusters
%    n_isolated        number of points not in clusters
%    n_pts             number of points per cluster
%    numclust(1,2,3)   number of singlet, double, multiple clusters, where
%                      singlet clusters include isolated points (see
%                      SRcluster.m for an equivalent definition)
%    singlet_faction   numclust(1) / sum(numclust)
%    sigma_actual      actual (computed) sigma of each cluster, that is, the
%                      standard deviation of the intracluster distances
%    indices_hull      cell array of boundary hull indices relative to XY per
%                      cluster
%    indConvex         as indices_hull but for convex hull boundaries
%    indCompact        as indices_hull but for the most compact boundaries
%    clust_width       maximum distance between vertices in each cluster
%    areas             area of each cluster
%    areasConvex       area of each cluster using the convex hull boundary
%    areasCompact      area of each cluster using the most compact boundary
%    equiv_radii       equivalent radius of each cluster
%    n_pts_per_area    number of points per area for clusters containing 3 or
%                      more points
%    perimeters        perimeter of each cluster
%    perimConvex       perimeter of each cluster using the convex hull boundary
%    perimCompact      perimeter of each cluster using the most compact
%                      boundary
%    compactness       4 pi area / perimeter^2 of each cluster
%                         [1 for a circle; 0 for a line segment]
%    circularity       4 pi areaCompact / perimConvex^2 of each cluster
%                         [1 for a circle; < 1 for non-circular, but relatively
%                          insensitive to irregular boundaries]
%    convexity         perimConvex / perimCompact of each cluster
%                         [1 for a convex object; < 1 if non-convex, e.g., one
%                          having an irregular boundary]
%    solidity          areaCompact / areaConvex of each cluster
%                         [1 for a solid object; < 1 if irregular boundary or
%                          containing holes]
%    min_c2c_dists     minimum center-to-center distances for each cluster
%                      with respect to all the others
%    min_e2e_dists     minimum edge-to-edge distances for each cluster convex
%                      hull with respect to all the others
%    min_c2c_dist      min(min_c2c_dists)
%    min_e2e_dist      min(min_e2e_dists)
%    nn_within_clust   nearest neighbor distances between points within
%                      clusters only
%
% Reference:
%    http://www.cyto.purdue.edu/cdroms/micro2/content/education/wirth10.pdf

% Created by
%    Michael J. Wester (2021)

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

   dim = size(XY, 2);

   nC = length(C);
   results.nC = nC;
   results.C  = C;
   results.XY = XY;

   xy_hull       = cell(1, nC);
   min_e2e_dists = zeros(1, nC);
   results.n_points     = size(XY, 1);
   results.n_clustered  = 0;
   results.n_isolated   = 0;
   results.n_pts        = zeros(1, nC);
   results.numclust     = zeros(1, 3);
   results.sigma_actual = zeros(1, nC);
   results.indices_hull = cell(1, nC);
   results.indConvex    = cell(1, nC);
   results.indCompact   = cell(1, nC);
   results.clust_width  = zeros(1, nC);
   results.areas        = zeros(1, nC);
   results.areasConvex  = zeros(1, nC);
   results.areasCompact = zeros(1, nC);
   results.equiv_radii  = zeros(1, nC);
   results.n_pts_per_area = [];
   results.perimeters   = [];
   results.perimConvex  = [];
   results.perimCompact = [];
   results.compactness  = [];
   results.circularity  = [];
   results.convexity    = [];
   results.solidity     = [];
   results.nn_within_clust = [];
   for i = 1 : nC
      xy = double(XY(C{i}, :));
      n_pts = size(xy, 1);
      results.n_pts(i) = n_pts;
      results.n_clustered = results.n_clustered + n_pts;
      if n_pts == 1
         results.numclust(1) = results.numclust(1) + 1;
         results.sigma_actual(i) = 0;
         results.indices_hull{i} = 1;
         results.indConvex{i}    = 1;
         results.indCompact{i}   = 1;
         results.clust_width(i)  = 0;
         results.areas(i)        = 0;
         results.areasConvex(i)  = 0;
         results.areasCompact(i) = 0;
         results.equiv_radii(i)  = 0;
         results.compactness(i)  = 1;
         results.circularity(i)  = 1;
         results.convexity(i)    = 1;
         results.solidity(i)     = 1;
      elseif n_pts == 2
         results.numclust(2) = results.numclust(2) + 1;
         results.sigma_actual(i) = std(pdist(xy));
         results.indices_hull{i} = [1, 2];
         results.indConvex{i}    = [1, 2];
         results.indCompact{i}   = [1, 2];
         results.areas(i)        = pdist(xy);
         results.areas(i)        = 0;
         results.areasConvex(i)  = 0;
         results.areasCompact(i) = 0;
         results.equiv_radii(i)  = pdist(xy) / 2;
         results.compactness(i)  = 0;
         results.circularity(i)  = 0;
         results.convexity(i)    = 1;
         results.solidity(i)     = 1;
         results.nn_within_clust = [results.nn_within_clust, ...
                                    repmat(pdist(xy), 1, 2)];
      else % n_pts >= 3
         results.numclust(3) = results.numclust(3) + 1;
         % The line below can sometimes cause MATLAB to crash inelegantly if
         % n_pts is large and may have some special size.
         %results.sigma_actual(i) = std(pdist(xy));
         % DoSigmaActual is a flag to allow this calculation to be performed,
         % noting the potential memory issues.
         if obj.DoSigmaActual
            results.sigma_actual(i) = std(pdist(xy));
         else
            results.sigma_actual(i) = 0;
         end
         try
            %[k, A] = convhull(xy(:, 1), xy(:, 2));
            % ShrinkFactor in boundary [DEFAULT = 0.5]
            %    = 0 gives the convex hull
            %    = 1 gives a compact boundary that envelops the points
            if dim == 2
               [k, A] = boundary(xy(:, 1), xy(:, 2), obj.ShrinkFactor);
               [kConvex,  AConvex]  = boundary(xy(:, 1), xy(:, 2), 0);
               [kCompact, ACompact] = boundary(xy(:, 1), xy(:, 2), 1);
            else
               [k, A] = boundary(xy(:, 1), xy(:, 2), xy(:, 3), ...
                                 obj.ShrinkFactor);
               [kConvex,  AConvex]  = ...
                  boundary(xy(:, 1), xy(:, 2), xy(:, 3), 0);
               [kCompact, ACompact] = ...
                  boundary(xy(:, 1), xy(:, 2), xy(:, 3), 1);
            end
            if isempty(k)
               k = 1;
               kConvex  = 1;
               kCompact = 1;
            end
         catch
            fprintf('boundary collinear (n_points = %d)\n', n_pts);
            %xy
            k = 1;
            kConvex  = 1;
            kCompact = 1;
            A = 0;
            AConvex  = 0;
            ACompact = 0;
         end
         results.indices_hull{i} = k;
         results.indConvex{i}    = kConvex;
         results.indCompact{i}   = kCompact;
         results.clust_width(i)  = smi_cluster.Clustering.vertex2vertex( ...
            xy(results.indices_hull{i}, :));
         results.areas(i)        = A;
         results.areasConvex(i)  = AConvex;
         results.areasCompact(i) = ACompact;
         if dim == 2
            % A = pi r^2
            results.equiv_radii(i)  = sqrt(A / pi);
         else
            % V (A) = 4/3 pi r^3
            results.equiv_radii(i)  = (3/4*A / pi)^(1/3);
         end
         if A > 0
            results.n_pts_per_area  = [results.n_pts_per_area, n_pts / A];
            if dim == 2
               perim = 0;
               perimConvex  = 0;
               perimCompact = 0;
               for j = 1 : length(k) - 1
                  perim = perim + pdist([xy(k(j), :); xy(k(j + 1), :)]);
               end
               for j = 1 : length(kConvex) - 1
                  perimConvex = perimConvex + ...
                     pdist([xy(kConvex(j), :); xy(kConvex(j + 1), :)]);
               end
               for j = 1 : length(kCompact) - 1
                  perimCompact = perimCompact + ...
                     pdist([xy(kCompact(j), :); xy(kCompact(j + 1), :)]);
               end
               results.perimeters      = [results.perimeters, perim];
               results.perimConvex     = [results.perimConvex,  perimConvex];
               results.perimCompact    = [results.perimCompact, perimCompact];
               results.compactness     = [results.compactness, ...
                                          4*pi*A / perim^2];
               results.circularity     = [results.circularity, ...
                                          4*pi*ACompact / perimConvex^2];
               results.convexity       = [results.convexity,   ...
                                          perimConvex / perimCompact];
               results.solidity        = [results.solidity ,   ...
                                          ACompact / AConvex];
            end
         end
      end
      xy_hull{i} = xy(results.indices_hull{i}, :);
      results.nn_within_clust = [results.nn_within_clust, ...
                                 smi_cluster.Clustering.nn_distances(xy)];
   end
   results.n_isolated = results.n_points - results.n_clustered;
   results.numclust(1) = results.numclust(1) + results.n_isolated;
   results.singlet_fraction = results.numclust(1) / sum(results.numclust);

   % The below can be an expensive operation timewise.
   if nC > 1 
      %results.min_c2c_dists = ...
      %   min(squareform(pdist(centers')) + 1.0e+10 * eye(nC));
      results.min_c2c_dists = smi_cluster.Clustering.nn_distances(centers');
      results.min_c2c_dist = min(results.min_c2c_dists);

      results.min_e2e_dists = zeros(1, nC);
      min_e2e_dist = 1.0e+10;
      for i = 1 : nC
         min_e2e_dists(i) = 1.0e+10;
         for j = 1 : i - 1
            e2e = smi_cluster.Clustering.edge2edge(xy_hull{i}, xy_hull{j});
            min_e2e_dist = min(min_e2e_dist, e2e);
            min_e2e_dists(i) = min(min_e2e_dists(i), e2e);
         end
         for j = i + 1 : nC
            e2e = smi_cluster.Clustering.edge2edge(xy_hull{i}, xy_hull{j});
            min_e2e_dist = min(min_e2e_dist, e2e);
            min_e2e_dists(i) = min(min_e2e_dists(i), e2e);
         end
         results.min_e2e_dists(i) = min_e2e_dists(i);
      end
      results.min_e2e_dist = min_e2e_dist;
   else
      results.min_c2c_dists = [];
      results.min_c2c_dist  = [];
      results.min_e2e_dists = [];
      results.min_e2e_dist  = [];
   end

end
