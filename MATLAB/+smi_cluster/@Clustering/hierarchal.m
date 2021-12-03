function [C, ptsI] = hierarchal(XY, E, minPts)
% Form clusters such that any point in a cluster is within E of some other
% point in the same cluster.
%
% Inputs:
%    XY       matrix of coordinates [N x 2]
%    E        cutoff distance
%    minPts   minimum number of points required for a cluster
%
% Outputs:
%    C        cell array of indices (wrt XY) per cluster [N x 1]
%    ptsI     indices (wrt XY) of isolated points

   if size(XY, 1) == 1
      C{1} = [1];
      ptsI = [];
      return;
   end

   Z = linkage(XY, 'single');
   %figure; dendrogram(Z);
   T = cluster(Z, 'Cutoff', E, 'Criterion', 'distance', 'Depth', 2);
   nC = max(T);

   % Remove clusters of size < minPts, taking the points to be isolated.
   C = [];
   j = 0;
   ptsI = [];
   for i = 1 : nC
      c = find(T == i);
      n = length(c);
      if n >= minPts
         j = j + 1;
         C{j} = c';
      else
         ptsI = [ptsI, c'];
      end
   end
   nC = j;

   ptsI = sort(ptsI);

%  % Check if the isolated points really are isolated or if they can be added
%  % to an existing cluster.  Do a single check here, although this might be
%  % iterated in the general case (or the algorithm rewritten entirely).
%  E2 = E^2;
%  not_isolated = [];
%  for i = 1 : length(ptsI)
%     I = ptsI(i);
%     p = XY(I, :);
%     MIN_d2 = 1.0e10;  
%     for j = 1 : nC
%        c = XY(C{j}, :);
%        min_d2 = min((c(:, 1) - p(1)).^2 + (c(:, 2) - p(2)).^2);
%        if min_d2 < MIN_d2
%           MIN_d2 = min_d2;
%           indx = j;
%        end
%     end
%     % Point is not isolated after all!  Add it into the nearest cluster.
%     if MIN_d2 < E2
%        C{indx} = [C{indx}, I];
%        not_isolated = [not_isolated, I];
%     end
%  end
%  % Removed non-isolated points (as determined above) from the list of
%  % isolated points.
%  ptsI = setdiff(ptsI, not_isolated);

end
