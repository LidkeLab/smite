function [N, M] = zNoll2NM(NollIndex)
% Convert Zernike index l into (n, m) (Noll ordering).

% Created by
%    Michael Wester, 2017, Lidkelab.

   N = ceil((-3 + sqrt(1 + 8*NollIndex)) / 2);
   M = NollIndex - N * (N + 1) / 2 - 1;
   if mod(N, 2) ~= mod(M, 2)
      M = M + 1;
   end
   if mod(NollIndex, 2) == 1
      M = -M;
   end

end
