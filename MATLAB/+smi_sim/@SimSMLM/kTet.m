function [x, y] = kTet(k, center, radius)
%Ktet produces a circle of k equally spaced points starting at a random place
% on the circumference.  Each segment centered around a point will subtend an
% angle of 2 pi/k.  The units of the output will be consistent with those of
% the inputs.
%
% INPUTS:
%    k        order of the k-tet (= number of points generated)
%    center   (x, y) coordinates of the circle's center
%    radius   radius of the circle
%
% OUTPUT:
%    x, y     coordinates of the localization computed

% Created by
%    Michael J. Wester (Lidkelab 2020)

   if k <= 0
      error('Ktet: k must be an integer > 0');
   end

   KtetAngle = (2 * pi) / k;
   % Start the k-tet at a random location along the circle's circumference.
   theta = 2 * pi * rand;
   x(1) = center(1) + radius * cos(theta);
   y(1) = center(2) + radius * sin(theta);

   for i = 2 : k
      theta = theta + KtetAngle;
      x(i) = center(1) + radius * cos(theta);
      y(i) = center(2) + radius * sin(theta);
   end

end
