function SMD_True = kDisc(k, center, radius)
%kDisc produces a disc of k points.
% The units of the output will be consistent with those of the inputs.
%
% INPUTS:
%    k            order of the k-disc (= number of points generated)
%    center       (x, y) coordinates of the circle's center [1 x 2] (pixel)
%    radius       radius of the circle (pixel)
%
% OUTPUT:
%    SMD_True     SMD structure containing:
%       X, Y         coordinates of the localizations computed

% Created by
%    Michael J. Wester (Lidkelab 2025)

   if k <= 0
      error('kDisc: k must be an integer > 0');
   end

   x = zeros(k, 1);
   y = zeros(k, 1);

   for i = 1 : k
      r = radius * rand;
      theta = 2 * pi * rand;
      x(i) = center(1) + r * cos(theta);
      y(i) = center(2) + r * sin(theta);
   end

   SMD_True = smi_core.SingleMoleculeData.createSMD();
   SMD_True.X = x;
   SMD_True.Y = y;
   SMD_True.Z = [];

end
