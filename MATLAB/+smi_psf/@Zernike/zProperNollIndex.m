function l_max = zProperNollIndex(l_max)
%zProperNollIndex Make sure that l_max includes both azimuthal terms (cos + sin)
% for a particular n and m.  For m ~= 0, the m corresponding to l_max and
% l_max - 1 (or l_max + 1) should agree in absolute value if both azimuthal
% terms are included (m == 0 corresponds to the case of no azimuthal terms).

% Created by
%    Michael Wester, 2017, Lidkelab.

   [n, m] = SMA_PSF.zernikeNoll2NM(l_max);
   if m ~= 0
      [nn, mm] = SMA_PSF.zernikeNoll2NM(l_max - 1);
      if abs(m) ~= abs(mm)
         l_max = l_max + 1;
      end
   end

end
