function [n, m] = zWyant2NM(l)
% Convert Zernike index l into (n, m) (Wyant ordering).

% Created by
%    Michael Wester, 2017, Lidkelab.

   n = floor(sqrt(l));
   mm = ceil((2*n - (l - n^2)) / 2);
   if mod(l - n^2, 2) == 0
      m =  mm;
   else
      m = -mm;
   end

end
