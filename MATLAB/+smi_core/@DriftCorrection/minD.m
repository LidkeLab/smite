%function [sumNND, X]  = minD_intra(Theta, X, T, Ndims, L_intra)
%function [sumNND, X2] = minD_inter(Theta, NS, X2, Ndims, L_inter)
function [sumNND, X] = minD(Theta, X, T, Ndims, L, NS)
% Sum of nearest neighbor distances for intra/inter-dataset drift correction.
%
% INPUTS:
%   Theta   approximating polynomial coefficients
%   X       localization coordinates (NX x Ndims) [pixel]
%   T       localization frame numbers (NX x Ndims)
%   Ndims   coordinate dimension (2 or 3)
%   L       threshold [pixel]
%   Y       localization coordinates for a 2nd dataset (NY x Ndims) [pixel]
%   NS      nearest neighbor searcher object referencing a 2nd static dataset
%           that will compared to X many times (Y); created by createns
%
% OUTPUTS:
%   sumNND  sum of nearest neighbor distances used as a cost estimate for the
%           optimizer
%   X       drift corrected coordinates (NX x Ndims) [pixel]

   if ~exist('NS', 'var') || isempty(NS)
      % --- Intra-dataset ---
      N = numel(Theta);

      m = N/Ndims;
      lo = 1;   hi = m;
      for j = 1:Ndims
         PX = Theta(lo : hi);
         % Assume the constant term of the polynomial is zero for the frames
         % within a dataset.  (This assumption considerably speeds up the
         % intra-dataset minimization.)
         X(:, j) = X(:, j) - polyval([PX', 0], T);

         lo = lo + m;
         hi = hi + m;
      end

      % D(:, 1) will be identically zero since here we are finding the nearest
      % neighbors of X with respect to itself.
      [~, D] = knnsearch(X, X, 'K', 2);
      D = D(:, 2);
   else
      % --- Inter-dataset ---
      % Theta(1) = DriftX, Theta(2) = DriftY, Theta(3) = DriftZ
      for j = 1:Ndims
         X(:, j) = X(:, j) - Theta(j);
      end

      [~, D] = knnsearch(NS, X);
      %[~, D] = knnsearch(X, Y);
   end

   % Make distances > L equal to L in the cost function to
   % de-emphasize widely separated nearest neighbors.  This can be be important
   % in sparse datasets where long distances may separate localizations.  The
   % nearest neighbor of a localization in one frame of a dataset may blink
   % off, resulting in a totally different nearest neighbor for the same
   % localization in the next frame.
   sumNND = sum(min(D, L));

   % In rare cases, D may contain NaNs, so here replace them with an average
   % value.  n is the number of NaNs found and the NaNs are removed from D
   % before recomputing sumNND.
   if isnan(sumNND)
      nans = isnan(D);
      D(nans) = [];
      n = sum(nans);
      % Note that sum(D(D <= L)) + L*sum(D > L) = sum(min(D, L)), then
      % sumNND = sum(min(D, L)) + n*mean(min(D, L))
      %        = sum(min(D, L)) + n * sum(min(D, L))/numel(min(D, L))
      %        = sum(min(D, L)) * (1 + n/numel(D));
      sumNND = sum(min(D, L)) * (1 + n/numel(D));
   end

end
