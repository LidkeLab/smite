function simkDiscs(obj, kk, radius_kDisc)
%simkDiscs generates 2D k-discs in the simulation region (units are pixels).
% See kDisc for further information.  The number of k-discs generated is based
% on the class property Rho (fluorophore density).
%
% INPUTS:
%    kk           order of the k-Discs
%    radius_kDisc radius of the circles (pixel)
%
% OUTPUT:
%    obj.SMD_True    SMD structure containing:
%       X, Y         coordinates of the localizations computed (pixel)

% Created by
%    Michael J. Wester (Lidkelab 2025)

   kDiscDensity = obj.Rho / kk;   % obj.Rho is in fluorophores/pixel
   n_kDiscs = round(kDiscDensity * obj.SZ^2);
   if obj.Verbose >= 1
      fprintf('Generating %d %d-Discs ...\n', n_kDiscs, kk);
   end
   center_kDisc = rand(1, 2) * obj.SZ;
   obj.SMD_True = obj.kDisc(kk, center_kDisc, radius_kDisc);
   for i = 2 : n_kDiscs
      center_kDisc = rand(1, 2) * obj.SZ;
      SMD_True_tmp = obj.kDisc(kk, center_kDisc, radius_kDisc);
      obj.SMD_True = ...
         smi_core.SingleMoleculeData.catSMD(obj.SMD_True, SMD_True_tmp, false);
   end
   obj.SMD_True.NDims = 2;
   obj.SMD_True.XSize = obj.SZ;
   obj.SMD_True.YSize = obj.SZ;

   % Apply labeling efficiency and generate blinks.
   obj.genModel();

end
