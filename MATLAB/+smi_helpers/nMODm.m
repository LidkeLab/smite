function r = nMODm(n, m)
% Modulus such that r is in [1, m] rather than [0, m - 1].
%
% Input:
%    n   integer
%    m   integer modulus
% Output:
%    r   remainder

   r = mod(n, m);
   if r <= 0
      r = r + m;
   end

end
