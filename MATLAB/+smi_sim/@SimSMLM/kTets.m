function SMD_True = kTets(obj, kk, radius_kTet)
%kTet generates 2D k-tets in the simulation region (units are pixels).
% See kTet for further information.  The number of k-tets generated is based on
% the class property Rho (fluorophore density).
%
% INPUTS:
%    kk           order of the k-tets
%    radius_kTet  radius of the circles
%
% OUTPUT:
%    SMD          SMD structure containing:
%       X, Y         coordinates of the localizations computed

% Created by
%    Michael J. Wester (Lidkelab 2021)

   kTetDensity = obj.Rho / kk;   % obj.Rho is in fluorophores/pixel
   n_kTets = round(kTetDensity * obj.SZ^2);
   if obj.Verbose >= 1
      fprintf('Generating %d %d-tets ...\n', n_kTets, kk);
   end
   center_kTet = rand(1, 2) * obj.SZ;
   SMD_True = obj.kTet(kk, center_kTet, radius_kTet);
   for i = 2 : n_kTets
      center_kTet = rand(1, 2) * obj.SZ;
      SMD_True_tmp = obj.kTet(kk, center_kTet, radius_kTet);
      SMD_True = smi_core.SingleMoleculeData.catSMD(SMD_True, SMD_True_tmp);
   end
   SMD_True.NDims = 2;
   SMD_True.XSize = obj.SZ;
   SMD_True.YSize = obj.SZ;

end
