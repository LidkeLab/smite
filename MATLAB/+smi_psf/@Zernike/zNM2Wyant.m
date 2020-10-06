function l = zNM2Wyant(n, m)
%zNM2Wyant Convert Zernike indices (n, m) into l (Wyant ordering).

% Created by
%    Michael Wester, 2017, Lidkelab.

   if m >= 0
      l = n^2 + 2*(n - m);
   else
      l = n^2 + 2*(n + m) + 1;
   end

end
