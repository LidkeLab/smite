function H = hopkinstat3(P,A,B,C,mm)
% Compute a Hopkin's statistic for the 3D particles P.
% Written by Michael Wester and Stanly Steinberg in 2008.
% The particle region is [0,A] x [0,B] x [0,C].
% m is the number of tests points.
% n is the number of particles.
%n = length(P);
n = numel(P) / 3;
m = min(mm, n);
% Create the indices of m test particles.
index = [0,m+1];
while ( (length(index) < m) | (index(1) < 1) | (index(end) > n) )
   index = unique(round(1/2+n*rand(m,1)));
end
T = P(index,:);
% Create m test points
S = rand(m,3)*diag([A,B,C]);
% Compute the minimum distance.
U = ones(1,m)*sqrt(A^2+B^2+C^2);
% W = ones(1,m)*sqrt(A^2+B^2+C^2);
W = U;
for k = 1:m
   for i = 1:n
      dist = sqrt((S(k,1)-P(i,1))^2 + (S(k,2)-P(i,2))^2 + (S(k,3)-P(i,3))^2);
      if dist > 0 
         U(k) = min(U(k), dist);
      end

      dist = sqrt((T(k,1)-P(i,1))^2 + (T(k,2)-P(i,2))^2 + (T(k,3)-P(i,3))^2);
      if dist > 0 
         W(k) = min(W(k), dist);
      end
   end
end
H = (sum(U.^3)/(sum(U.^3)+sum(W.^3)));

end
