function plot_voronoi3(X, Y, Z, v, c, rho, str, dense, ptIDs)
%plot_voronoi3 plots the Voronoi diagram corresponding to (X, Y, Z), coloring
% the cells according to the density rho.
%
% INPUTS:
%    X, Y, Z point coordinates, each [N x 1]
%    v       vertices [N x 2]
%    c       vertex indices cooresponding to each cell [Nc x 1]
%    rho     density of each cell [Nc x 1]
%    str     figure title
%    dense   indices of those cells >= prescribed density criterion
%    ptIDs   [OPTIONAL] write out point IDs if true [default = false]

% Created by
%    Michael J. Wester (2019)

   n_XY = length(X);
   n_v  = size(v, 1);
   n_c  = length(c);

   delta = 0.001 * (max(X) - min(X));
   figure();
   hold on
   dt = delaunayTriangulation(X, Y, Z);
   tetramesh(dt);
   limits = axis;
   colormap(jet);
   caxis([min(rho), max(rho)]);
   for i = 1 : n_c
      c_i = c{i};
      if all(c_i ~= 1)
         fill3(v(c_i, 1), v(c_i, 2), v(c_i, 3), rho(i));
         %fill(v(c_i, 1), v(c_i, 2), floor((i - 1)/(n_c - 1) * 255 + 1));
      end
   end
   plot3(X, Y, Z, 'k.', 'MarkerSize', 10);
   plot3(X(dense), Y(dense), Z(dense), 'r.', 'MarkerSize', 10);
   if exist('ptIDs', 'var') & ptIDs
      for i = 1 : n_XY
         text(X(i) + delta, Y(i), Z(i), sprintf('%d', i), 'Color', 'k');
      end
   end
%  plot(v(:, 1), v(:, 2), 'g.', 'MarkerSize', 10);
%  for i = 2 : n_v
%     text(v(i, 1) + delta, v(i, 2), sprintf('%d', i), 'Color', 'g');
%  end
   axis(limits);
   cb = colorbar;
   cb.Label.String = 'relative density';
   title(str);
   hold off

end
