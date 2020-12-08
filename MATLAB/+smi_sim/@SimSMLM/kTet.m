function SMD_True = kTet(k, center, radius, startAngle)
%ktet produces a circle of k equally spaced points starting at a random place
% on the circumference unless the optional argument startAngle is provided.
% Each segment centered around a point will subtend an angle of 2 pi/k.  The
% units of the output will be consistent with those of the inputs.
%
% INPUTS:
%    k            order of the k-tet (= number of points generated)
%    center       (x, y) coordinates of the circle's center [1 x 2]
%    radius       radius of the circle
%    startAngle   [OPTIONAL] starting angle for the first localization (radian)
%
% OUTPUT:
%    SMD          SMD structure containing:
%       X, Y         coordinates of the localizations computed

% Created by
%    Michael J. Wester (Lidkelab 2020)

   if k <= 0
      error('kTet: k must be an integer > 0');
   end

   x = zeros(k, 1);
   y = zeros(k, 1);

   kTetAngle = (2 * pi) / k;
   % Start the k-tet at a random location along the circle's circumference.
   if exist('startAngle', 'var')
      theta = startAngle;
   else
      theta = 2 * pi * rand;
   end
   x(1) = center(1) + radius * cos(theta);
   y(1) = center(2) + radius * sin(theta);

   for i = 2 : k
      theta = theta + kTetAngle;
      x(i) = center(1) + radius * cos(theta);
      y(i) = center(2) + radius * sin(theta);
   end

   SMD_True = smi_core.SingleMoleculeData.createSMD();
   SMD_True.X = x;
   SMD_True.Y = y;
   SMD_True.Z = [];

end
