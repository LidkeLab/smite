function simkTets(obj, kk, radius_kTet)
%simkTet generates 2D k-tets in the simulation region (units are pixels).
% See kTet for further information.  The number of k-tets generated is based on
% the class property Rho (fluorophore density).
%
% INPUTS:
%    kk           order of the k-tets
%    radius_kTet  radius of the circles (pixel)
%
% OUTPUT:
%    obj.SMD_True    SMD structure containing:
%       X, Y         coordinates of the localizations computed (pixel)

% Created by
%    Michael J. Wester (Lidkelab 2021)

   kTetDensity = obj.Rho / kk;   % obj.Rho is in fluorophores/pixel
   n_kTets = round(kTetDensity * obj.SZ^2);
   if obj.Verbose >= 1
      fprintf('Generating %d %d-tets ...\n', n_kTets, kk);
   end
   center_kTet = rand(1, 2) * obj.SZ;
   obj.SMD_True = obj.kTet(kk, center_kTet, radius_kTet);
   for i = 2 : n_kTets
      center_kTet = rand(1, 2) * obj.SZ;
      SMD_True_tmp = obj.kTet(kk, center_kTet, radius_kTet);
      obj.SMD_True = ...
         smi_core.SingleMoleculeData.catSMD(obj.SMD_True, SMD_True_tmp);
   end
   obj.SMD_True.NDims = 2;
   obj.SMD_True.XSize = obj.SZ;
   obj.SMD_True.YSize = obj.SZ;

end
