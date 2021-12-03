function [nC, C] = cluster_HSET(obj, SMD, minPts)
% Perform the clustering implied by H-SET, in which the nodes that would be
% combined in normal H-SET are taken to be clusters here.
%
% INPUTS:
%    obj         various properties used by the algorithms
%       PixelSize    linear dimension of a pixel (nm)
%       Sigma_Reg    pre-computed image registration error (nm)
%    SMD         SMD structure with data for H-SET clustering
%       X, Y         X and Y localization coordinates     (pixel)
%       X_SE, Y_SE   X and Y localization standard errors (pixel)
%    minPts      minimum number of points allowed in a cluster
%
% OUTPUTS:
%    nC          number of clusters found
%    C           cell array of XY indices forming each cluster [nC x 1]

% Created by
%    Michael J. Wester (2021)

   XY    = double([ SMD.X, SMD.Y ]) .* obj.PixelSize;
   sigma = double([ SMD.X_SE, SMD.Y_SE ]) .* obj.PixelSize;
   [XY_new, sigma_new, combined, SRsave] = ...
      obj.clusterSR(XY, sigma, obj.Sigma_Reg);

   nC = 0;
   C = {};
   for i = 1 : numel(combined)
      if numel(combined{i}) >= minPts
         nC = nC + 1;
         C{nC} = combined{i};
      end
   end
   
end
