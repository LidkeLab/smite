function min_d = nn_distances(xy)
%nn_distances is the minimum nearest neighbor distances from each point in xy
% to the other points.
%
% INPUTS:
%    xy      n x 2 (or n x 3) set of coordinates.
%
% OUTPUTS:
%    min_d   minimum nearest neighbor distance for each point in xy

% Created by
%    Michael J. Wester (2019)

   %x_r = xy(:, 1);
   %y_r = xy(:, 2);
   %min_d = min(squareform(pdist([x_r, y_r])) + 1.0e+10 * eye(numel(x_r)));
   % Below is much faster, less memory intensive and works in 3D as well as 2D!
   min_d = [];
   if size(xy, 1) > 1
      [~, min_d] = knnsearch(xy, xy, 'K', 2);
      min_d = min_d(:, 2)';
   end

end
