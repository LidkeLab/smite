function [XY_combined, sigma_combined, nodes_combined] = ...
   hierarchalSingleLabel(XY, sigma, Sigma_Reg, LoS)
% Combine multiple clustered points into single labels when appropriate via a
% top-down descent through a hierarchal dendrogram relationship between points.
% n is the original number of points and m is the dimension.
% n' is the final number of points after combinations have occurred.
%
% Inputs:
%    XY          n x m matrix of coordinates (nm)
%    sigma       n x m matrix of position uncertainties (1 std deviation) (nm)
%    Sigma_Reg   1 x m array of registration error (1 standard deviation) (nm)
%    LoS         level of significance [0 <= LoS <= 1 with 0.01 typical]
%
% Outputs:
%    XY_combined      n' x m final coordinate matrix (nm)
%    sigma_combined   n' x m final position uncertainty matrix (nm)
%    nodes_combined   cell array of indices of combined points per cluster

   XY_combined    = XY;
   sigma_combined = sigma;
   nodes_combined = {};

   n_pts   = size(XY, 1);
   n_nodes = n_pts - 1;   % does not include leaf nodes

   if n_pts <= 1
      return
   end

   Z = linkage(XY, 'single');
   %figure(); dendrogram(Z);

   % Find all the leaf nodes contained by each composite node by parsing the
   % Z matrix (see MATLAB linkage documentation).  Note that composite nodes
   % are indexed as node # - n_pts.  E.g., if n_pts = 10, then Z(2, :) is
   % composite node 2 and overall node 12.  Z(2, 1:2) are the nodes (in
   % overall node numbering) that are joined by node 12 which, for example,
   % might be 7 (a leaf node since it is <= 10) and 11 (a composite node whose
   % components are given in Z(1, 1:2)).
   cn = cell(1, n_nodes);
   for i = 1 : n_nodes
      l1 = Z(i, 1);   % child node 1
      l2 = Z(i, 2);   % child node 2
      if l1 <= n_pts
         leaf_nodes = l1;
      else
         leaf_nodes = cn{l1 - n_pts};
      end
      if l2 <= n_pts
         leaf_nodes = [leaf_nodes, l2];
      else
         leaf_nodes = [leaf_nodes, cn{l2 - n_pts}];
      end
      cn{i} = sort(leaf_nodes);
   end

   k = 0;
   % Start with the top-level node and descend down through the tree (the tree
   % is taken to have its root at the top and its leaves at the bottom).
   for i = n_nodes : -1 : 1
      LN = cn{i};
      if ~isempty(LN)
         [Pvalue, XY_wm, sigma_wm] = ...
            smi_cluster.Clustering.singleLabelTest(XY(LN, :), sigma(LN, :), ...
                                                   Sigma_Reg);
         % Delete leaf nodes that are part of a composite single label.
         % Retain the first leaf node to hold the new information.
         if Pvalue > LoS
            cn{i} = LN(1);
            LN_deleted = LN(2 : end);
            % To be deleted at the end, but in order to keep the numbering
            % unchanged within the loop, set to a fantasy value for now.
            XY_combined(LN_deleted, :)    = NaN;
            sigma_combined(LN_deleted, :) = NaN;
            % Replace XY and sigma for this first leaf node with the weighted
            % means of the collapsed cluster.
            XY_combined(LN(1), :)    = XY_wm;
            sigma_combined(LN(1), :) = sigma_wm;
            k = k + 1;
            nodes_combined{k} = LN;
            for j = 1 : i - 1
               % cn{j} = setdiff(cn{j}, LN);
               % --->
               % cn{j} = cn{j}(~ismember(cn{j}, LN));
               % Optimized version of the line above: 
               cn{j} = ...
cn{j}(~smi_cluster.Clustering.my_ismemberBuiltinTypes(cn{j}, LN));
            end
         end
      end
   end

   XY_combined(isnan(XY_combined(:, 1)), :) = [];
   sigma_combined(isnan(sigma_combined(:, 1)), :) = [];

end
