function SMD = inflateSE(SMD, SEAdjust)
%inflateSE: Inflate standard errors.
%
% INPUTS:
%    SMD        Single Molecule Data structure
%    SEAdjust   standard error inflation amount (pixel) [DEFAULT = 0]
%
% OUTPUT:
%    SMD        modified Single Molecule Data structure

% Created by
%    David J. Schodt and Michael J. Wester (5/24/2022)

if ~exist('SEAdjust', 'var')
   SEAdjust = 0;
end

SMD.X_SE = SMD.X_SE + SEAdjust;
SMD.Y_SE = SMD.Y_SE + SEAdjust;

end
