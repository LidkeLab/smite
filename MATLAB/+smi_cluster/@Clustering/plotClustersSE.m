function clusterFig = plotClustersSE(obj, SMD, C, centers, ptsI, txt, ...
                                     options)
% Plot and label the 2D clusters.
%
% INPUTS:
%    SMD          (x, y) coordinates in the format
%                     (1) SMD structure: SMD.X, SMD.Y, SMD.X_SE, SMD.Y_SE (pixel)
%                     (2) N x 4 array of coordinates and precsions (nm):
%                         X, Y, X_SE, Y_SE
%    C            cell array of XY indices forming each cluster [nC x 1]
%    centers      coordinates of the center of each cluster [nC x n_dim]
%    ptsI         indices of points not found in any cluster
%    txt          descriptive text to add to the plot's title
%    options      'L': label each cluster
%                 'O': outline each cluster
%                 'P': print the size of each cluster
%                 '1': use only one color for all clusters
% OUTPUT:
%    clusterFig   figure handle

   % Check for SMD structure; in this situation, assume the units are pixels.
   SMDstruct = false;
   if isfield(SMD, 'X') & isfield(SMD, 'Y')
      SMDstruct = true;
      XY = [SMD.X, SMD.Y] .* obj.PixelSize;
      XY_SE = [SMD.X_SE, SMD.Y_SE] .* obj.PixelSize;
   % Check for N x 2 or N x 3 matrix of coordinates.
   elseif ismatrix(SMD) && size(SMD, 2) == 4
      XY = SMD(:, 1:2);
      XY_SE = SMD(:, 3:4);
   else
      error('Unrecognized coordinate input!');
   end

   if exist('options', 'var')
      labeling = strfind(options, 'L');
      outlines = strfind(options, 'O');
      printing = strfind(options, 'P');
      onecolor = strfind(options, '1');
   else
      labeling = false;
      outlines = true;
      printing = false;
      onecolor = false;
   end

   colors = 'rgbcmy';
   n_colors = length(colors);

   nC = length(C);
   n_isolated = length(ptsI);
   n_points = size(XY, 1);
   n_clustered = 0;

   clusterFig = figure('Visible', 'off');
   hold on
   fprintf('\n');
   theta = 0 : 0.05 : 2*pi + 0.05;
   for i = 1 : nC
      j = smi_helpers.nMODm(i, n_colors);
      if onecolor
         color = 'm';
      else
         color = colors(j);
      end
      n_pts = length(C{i});
      n_clustered = n_clustered + n_pts;
      if printing
         fprintf('cluster %d has %2d points (%s)\n', i, n_pts, color);
      end
      %plot(XY(C{i}, 1), XY(C{i}, 2), [color, '.'], 'MarkerSize', 12);
      plot(XY(C{i}, 1), XY(C{i}, 2), [color, '.']);
      for j = 1 : n_pts
         r = sqrt(2*(XY_SE(C{i}(j), 1)^2 + XY_SE(C{i}(j), 2)^2));
         x = XY(C{i}(j), 1) + r*cos(theta);
         y = XY(C{i}(j), 2) + r*sin(theta);
         plot(x, y, [color, '-']);
      end
      if labeling
         text(centers(1, i), centers(2, i), sprintf('%d', i));
      end
      if outlines
         if length(C{i}) <= 2
            plot(XY(C{i}, 1), XY(C{i}, 2), [color, '-'], 'LineWidth', 2);
         else %if length(C{i}) >= 3
            xy = double(XY(C{i}, :));
            %k = convhull(XY(C{i}, 1), XY(C{i}, 2));
            k = boundary(xy(:, 1), xy(:, 2), obj.ShrinkFactor);
            k = C{i}(k);
            plot(XY(k, 1), XY(k, 2), [color, '-'], 'LineWidth', 2);
         end
      end
   end
   if n_isolated > 0
      %plot(XY(ptsI, 1), XY(ptsI, 2), 'k.', 'MarkerSize', 12);
      plot(XY(ptsI, 1), XY(ptsI, 2), 'k.');
      for j = 1 : n_isolated
         r = sqrt((XY_SE(ptsI(j), 1)^2 + XY_SE(ptsI(j), 2)^2) / 2);
         x = XY(ptsI(j), 1) + r*cos(theta);
         y = XY(ptsI(j), 2) + r*sin(theta);
         plot(x, y, 'k-');
      end
   end
   title(sprintf('%s (points = %d, clusters = %d, clustered %% = %.3f)', ...
                 regexprep(txt, '_', '\\_'), n_points, nC,               ...
                 n_clustered/n_points * 100));
   hold off

end
