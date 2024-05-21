function max_d = vertex2vertex(hull)
%vertex2vertex is the maximum vertex-to-vertex distance within the hull.
%
% INPUTS:
%    hull           (x, y {, z}) coordinates of the hull defining the
%                   boundaries of a cluster of points [N x dim]
%
% OUTPUT:
%    max_d          maximum vertex-to-vertex distance within the hull

% Created by
%    Michael J. Wester (2024)

   dim = size(hull, 2);
   max_d2 = -1.0;
   n1 = size(hull, 1);
   for i = 1 : n1
      d2 = 0;
      for j = 1 : dim
         d2 = d2 + (hull(i, j) - hull(:, j)).^2;
      end
      max_d2 = max(max_d2, max(d2));
   end
   max_d = sqrt(max_d2);

end
