function min_d = edge2edge(hull1, hull2)
%edge2edge is the minimum edge-to-edge distance between hull 1 and hull 2.
%
% INPUTS:
%    hull1, hull2   (x, y {, z}) coordinates of the two hulls defining the
%                   boundaries of a cluster of points [N x dim]
%
% OUTPUT:
%    min_d          minimum edge-to-edge distance between hull 1 and hull 2

% Created by
%    Michael J. Wester (2019)

   dim = size(hull1, 2);
   min_d2 = 1.0e+10;
   n1 = size(hull1, 1);
   %n2 = size(hull2, 1);
   for i = 1 : n1
      %d2 = (hull1(i, 1) - hull2(:, 1)).^2 + (hull1(i, 2) - hull2(:, 2)).^2;
      d2 = 0;
      for j = 1 : dim
         d2 = d2 + (hull1(i, j) - hull2(:, j)).^2;
      end
      min_d2 = min(min_d2, min(d2));
   end
   min_d = sqrt(min_d2);

end
