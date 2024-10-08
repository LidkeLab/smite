function clusterFig = plotClusters3(obj, SMD, C, centers, ptsI, txt, options)
% Plot and label the 3D clusters.
%
% INPUTS:
%    SMD         (x, y) coordinates in the format
%                    (1) SMD structure: SMD.X, SMD.Y with optional SMD.X_SE,
%                        SMD.Y_SE (only needed for H-SET clustering) (pixel)
%                    (2) N x 2 or N x 3 array of coordinates (nm)
%    C            cell array of XY indices forming each cluster [nC x 1]
%    centers      coordinates of the center of each cluster [nC x n_dim]
%    ptsI         indices of points not found in any cluster
%    txt          descriptive text to add to the plot's title
%    options      [default = 'O'] string containing some of the following:
%                 'L': label each cluster
%                 'O': outline each cluster
%                 'P': print the size of each cluster
%                 '1': use only one color for all clusters
% OUTPUT:
%    clusterFig   figure handle

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
   if printing
      fprintf('\n');
   end
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
      plot3(XY(C{i}, 1), XY(C{i}, 2), XY(C{i}, 3), [color, '.'], ...
            'MarkerSize', 12);
      if labeling
         text(centers(1, i), centers(2, i), centers(3, i), sprintf('%d', i));
      end
      if outlines
         if length(C{i}) <= 3
            plot3(XY(C{i}, 1), XY(C{i}, 2), XY(C{i}, 3), [color, '-'], ...
                  'LineWidth', 2);
         else %if length(C{i}) >= 4
            xy = double(XY(C{i}, :));
            %k = convhull(XY(C{i}, 1), XY(C{i}, 2), XY(C{i}, 3));
            k = boundary(xy(:, 1), xy(:, 2), xy(:, 3), obj.ShrinkFactor);
            k = C{i}(k);
            plot3(XY(k, 1), XY(k, 2), XY(k, 3), [color, '-'], 'LineWidth', 2);
            %trisurf(k, XY(:, 1), XY(:, 2), XY(:, 3), 'FaceColor', color);
         end
      end
   end
   if n_isolated > 0
      plot3(XY(ptsI, 1), XY(ptsI, 2), XY(ptsI, 3), 'k.', 'MarkerSize', 12);
   end
   title(sprintf('%s (clusters = %d, clustered %% = %.3f)', ...
                regexprep(txt, '_', '\\_'), nC, n_clustered/n_points * 100));
   hold off

end
