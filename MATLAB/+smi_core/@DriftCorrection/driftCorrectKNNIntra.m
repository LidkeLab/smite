function [SMD, Statistics] = driftCorrectKNNIntra(obj, SMD, iDataset)
%driftCorrectKNNIntra calculates intra-dataset drift directly from X,Y{,Z} coordinates
% by fitting a polynomial depending on time (i.e., frame number) to the frames
% with each dataset (intra-dataset), and fitting constant shifts between
% datasets (inter-dataset).  Fitting is done via performing fminsearch on the
% (weighted) sums of nearest neighbor distances.  Intra-dataset portion.
%
% INPUTS:
%   iDataset:    Dataset index
%   SMD:         A structure with fields:
%      X              x coordinates (Nx1) where N is total number of points
%      Y              y coordinates (Nx1)
%      Z              z coordinates (Nx1) [OPTIONAL]
%                  Note: X and Y are typically in pixels, Z in um
%      DatasetNum     dataset number from which localization originates (Nx1)
%      FrameNum       frame   number from which localization originates (Nx1)
%      NDatasets      number of datasets
%      NFrames        number of frames in each dataset
%   obj:           [class properties]
%                     optimization parameters with the following fields:
%      L_intra        intra-dataset threshold (Default = 1 pixel)
%      L_inter        inter-dataset threshold (Default = 2 pixels)
%      PixelSizeZUnit pixel size in um (units of Z; only needed for 3D)
%                     (Default = 0.1)
%      PDegree        degree of the intra-dataset drift correction fitting
%                     polynomial (Default = 1)
%      TolFun_intra   termination tolerance on the intra-dataset function value
%                     (Default = 1e-2)
%      TolX_intra     termination tolerance on the intra-dataset fitting
%                     polynomial (Default = 1e-4)
%      TolFun_inter   termination tolerance on the inter-dataset function value
%                     (Default = 1e-2)
%      TolX_inter     termination tolerance on the inter-dataset fitting
%                     polynomial (Default = 1e-4)
%      Init_inter     inter-dataset initialization with respect to the previous
%                     dataset; the value should be either 0 (no initial drift),
%                     1 (initial drift of the previous dataset) or SMD.NFrames
%                     (final drift); zero or initial drift should work well
%                     with brightfield registration, while final drift works
%                     well generally (but the optimization process may not
%                     converge quite as quickly) (Default = SMD.NFrames)
%
% OUTPUTS:
%   SMD:         SMD data structure with updated fields:
%      X              drift corrected x coordinates (Nx1)
%      Y              drift corrected y coordinates (Nx1)
%      Z              drift corrected z coordinates (Nx1) [OPTIONAL]
%      DriftX         found x drift (NFrames x NDatasets)
%      DriftY         found y drift (NFrames x NDatasets)
%      DriftZ         found z drift (NFrames x NDatasets) [OPTIONAL]
%   Statistics:  statistical information about the algorithm performance
%                including various input parameters above and ...:
%      NDatasets          internal number of datasets
%      NFrames            internal number of frames per dataset
%      Intra_iterations   intra-dataset number of fminsearch iterations
%      Intra_funcCount    intra-dataset number of fminsearch function evals
%      Intra_elapsedTime  intra-dataset elapsed time for drift correction
%      Inter_iterations   inter-dataset number of fminsearch iterations
%      Inter_funcCount    inter-dataset number of fminsearch function evals
%      Inter_elapsedTime  inter-dataset elapsed time for drift correction
%
% CITATION:
%    "Robust, Fiducial-Free Drift Correction for Super-resolution Imaging"

% Created by
%    Modelz... Bewerdorf's group
%    Farzin Farzam  (Lidke Lab 2017) [original version of driftCorrect2D]
%    Michael J. Wester and Keith Lidke (Lidke Lab 2017-2019) [knnsearch]
%       (completely rewrote the code)

   if isfield(SMD, 'Z') && numel(SMD.Z) == numel(SMD.X)
      Ndims = 3;
   else
      Ndims = 2;
   end

   DriftParams.L_intra        = obj.L_intra;
   DriftParams.L_inter        = obj.L_inter;
   DriftParams.PixelSizeZUnit = obj.PixelSizeZUnit;
   DriftParams.PDegree        = obj.PDegree;
   DriftParams.TolFun_intra   = obj.TolFun_intra;
   DriftParams.TolX_intra     = obj.TolX_intra;
   DriftParams.TolFun_inter   = obj.TolFun_inter;
   DriftParams.TolX_inter     = obj.TolX_inter;
   DriftParams.Init_inter     = obj.Init_inter;
   if ~isempty(obj.NDatasets)
      DriftParams.NDatasets   = obj.NDatasets;
   end
   if ~isempty(obj.NFrames)
      DriftParams.NFrames     = obj.NFrames;
   end

   % Initialize various parameters.
   % PixelSizeZUnit is needed for 3D to convert Z into the same units as X
   % and Y.  PixelSizeZUnit is the X/Y pixel size in units of um.
   PixelSizeZUnit = DriftParams.PixelSizeZUnit;
   % L_intra is the intra-dataset threshold used in the cost function
   % computation to limit spurious nearest neighbor distances.  L_intra is
   % in units of pixels.
   L_intra = DriftParams.L_intra;
   % L_inter is the inter-dataset threshold used in the cost function
   % computation to limit spurious nearest neighbor distances.  L_inter is
   % in units of pixels.
   L_inter = DriftParams.L_inter;
   PDegree      = DriftParams.PDegree;
   TolFun_intra = DriftParams.TolFun_intra;
   TolX_intra   = DriftParams.TolX_intra;
   TolFun_inter = DriftParams.TolFun_inter;
   TolX_inter   = DriftParams.TolX_inter;
   Init_inter   = DriftParams.Init_inter;

   Statistics.Ndims          = Ndims;
   Statistics.PixelSizeZUnit = PixelSizeZUnit;
   Statistics.NDatasets      = SMD.NDatasets;
   Statistics.NFrames        = SMD.NFrames;
   Statistics.L_intra        = L_intra;
   Statistics.L_inter        = L_inter;
   Statistics.PDegree        = PDegree;
   Statistics.Init_inter     = Init_inter;
   Statistics.TolFun_intra   = TolFun_intra;
   Statistics.TolX_intra     = TolX_intra;
   Statistics.TolFun_inter   = TolFun_inter;
   Statistics.TolX_inter     = TolX_inter;

   SMD.DriftX = zeros(SMD.NFrames, SMD.NDatasets, 'single');
   SMD.DriftY = zeros(SMD.NFrames, SMD.NDatasets, 'single');
   if Ndims == 3
      SMD.DriftZ = zeros(SMD.NFrames, SMD.NDatasets, 'single');
   end

   % ---------- Intra-Dataset drift correction --------------------------------

   tic;
   % These options will be used by fminsearch below.
   options = optimset('TolFun', TolFun_intra, 'TolX', TolX_intra);

   % Establish an indexing array to record rearrangements of points into
   % datasets, so that this reordering can be later undone.
   N = numel(SMD.X);
   idx = zeros(N, 1);
   base = 0;

   %SMRS = cell(1, SMD.NDatasets);
   % Count the number of iterations and function calls.
   it = 0;   fc = 0;
   % Note that variables like X, Y, Z are vectors, while variables like XY are
   % n x Ndims matrices.  The basic reason for combining vectors into matrices
   % is to allow easier generalization to 3D.
   for i = 1:SMD.NDatasets
      %mask = SMD.DatasetNum == i;
      mask = SMD.DatasetNum == iDataset;
      n = sum(mask);
      XY = zeros(n, Ndims, 'single');
      XY(:, 1) = SMD.X(mask);
      XY(:, 2) = SMD.Y(mask);
      if Ndims == 3
         % Assume z values are given in um.  Convert these into the same units
         % as the x and y values, which are assumed to be given in pixels.
         % The z values are converted back to um on output.
         XY(:, 3) = SMD.Z(mask) ./ PixelSizeZUnit;
      end
      FrameNum = single(SMD.FrameNum(mask));

      idx(mask) = base + (1 : n);
      base = base + n;

      % Initialize polynomial fitting function, then optimize over the dataset.
      % P0 and P are packed vectors containing coefficients for all Pdegree
      % interpolations over all dimensions as required by fminsearch.  Note
      % that the polynomial is of PDegree for the drift correction.
      P0 = zeros(Ndims*PDegree, 1);
      % Test for really sparse datasets.
      if n <= 1
         P = P0;
         XYC = XY;
      else
         [P, ~, exitflag, output] = ...
            fminsearch(@minD_intra, P0, options, XY, FrameNum, Ndims, L_intra);
         if exitflag ~= 1
            fprintf( ...
               'driftCorrectKNN fminsearch on minD_intra exitflag = %d\n', ...
                     exitflag);
         end
         it = it + output.iterations;
         fc = fc + output.funcCount;

         [~, XYC] = minD_intra(P, XY, FrameNum, Ndims, L_intra);
      end

      PX = P(1           : PDegree);
      PY = P(PDegree + 1 : 2*PDegree);

      % Values corrected for drift.
      %SMRS{i}.XY = XYC;
      %SMRS{i}.n  = n;
      SMRS{iDataset}.XY = XYC;
      SMRS{iDataset}.n  = n;

      range = double(1:SMD.NFrames);
      SMD.DriftX(:, i) = polyval([PX', 0], range);
      SMD.DriftY(:, i) = polyval([PY', 0], range);

      if Ndims == 3
         PZ = P(2*PDegree + 1 : end);
         SMD.DriftZ(:, i) = polyval([PZ', 0], range);
      end
   end
   Statistics.Intra_iterations  = it;
   Statistics.Intra_funcCount   = fc;
   Statistics.Intra_elapsedTime = toc;

   obj.idx = [obj.idx; idx + numel(obj.idx)];
   obj.SMRS{iDataset} = SMRS{iDataset};

   % ---------- Inter-Dataset drift correction --------------------------------

   % --------------------------------------------------------------------------

end

% =============================================================================

function [sumNND, X] = minD_intra(Theta, X, T, Ndims, L_intra)
% Sum of nearest neighbor distances for intra-dataset drift correction.

   N = numel(Theta);

   m = N/Ndims;
   lo = 1;   hi = m;
   for j = 1:Ndims
      PX = Theta(lo : hi);
      % Assume the constant term of the polynomial is zero for the frames
      % within a dataset.  (This assumption considerably speeds up the
      % intra-dataset minimization.)
      X(:, j) = X(:, j) + polyval([PX', 0], T);

      lo = lo + m;
      hi = hi + m;
   end

   % D(:, 1) will be identically zero since here we are finding the nearest
   % neighbors of X with respect to itself.
   [~, D] = knnsearch(X, X, 'K', 2);
   D2 = D(:, 2);

   % Make distances > L_intra equal to L_intra in the cost function to
   % de-emphasize widely separated nearest neighbors.  This can be be important
   % in sparse datasets where long distances may separate localizations.  The
   % nearest neighbor of a localization in one frame of a dataset may blink
   % off, resulting in a totally different nearest neighbor for the same
   % localization in the next frame.
   sumNND = sum(min(D2, L_intra));

   % In rare cases, D2 may contain NaNs, so here replace them with an average
   % value (modeled on the code in minD_inter).
   if isnan(sumNND)
      nans = isnan(D2);
      D2(nans) = [];
      n = sum(nans);
      sumNND = sum(min(D2, L_intra)) * (1 + n/numel(D2));
   end

end
