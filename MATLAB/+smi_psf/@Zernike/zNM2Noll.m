function NollIndex = zNM2Noll(N, M)
%zNM2Noll Convert Zernike indices (n, m) into l (Noll ordering).

% Created by
%    Michael Wester, 2017, Lidkelab.

   mm = abs(M);
   % l = n * (n + 1) / 2 + 1; if mm >= 2; l = l + mm - 1; end
   NollIndex = N * (N + 1) / 2 + 1 + max(0, mm - 1);
   if (M > 0 & mod(N, 4) >= 2) | (M < 0 & mod(N, 4) <= 1)
      NollIndex = NollIndex + 1;
   end

end
