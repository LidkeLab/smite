function [SMD, Statistics] = driftCorrectKNN(obj, SMD)
%driftCorrectKNN calculates the drift directly from X,Y{,Z} coordinates
% by fitting a polynomial depending on time (i.e., frame number) to the frames
% with each dataset (intra-dataset), and fitting constant shifts between
% datasets (inter-dataset).  Fitting is done via performing fminsearch on the
% (weighted) sums of nearest neighbor distances.
%
% INPUTS:
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
%      Verbose        verbosity level (Default = 1)
%      NDatasets      [OPTIONAL] override the collected value.  This causes the
%                     dataset/frame numbering to be reorganized internally as
%                     specified by the user
%      NFrames        [OPTIONAL] override the collected value.  See above
%      NOTES: Only one of NDatasets or NFrames needs to be specified.  These
%             numbers must evenly divide the total number of frames.  Better
%             results can sometimes occur by increasing the number of datasets
%             or decreasing the number of frames per dataset up to some limit
%             when the datasets become too sparse.  Init_inter will be
%             automatically changed from SMD.NFrames (if so specified) to
%             NFrames.
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
%      NDatasets_C        original (collected) number of datasets
%      NFrames_C          original (collected) number of frames per dataset
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

   NDatasets_C = SMD.NDatasets;
   NFrames_C   = SMD.NFrames;

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
   DriftParams.Verbose        = obj.Verbose;

   % Initialize various parameters, either provided by the user or defaults.
%  if exist('DriftParams', 'var')
      % PixelSizeZUnit is needed for 3D to convert Z into the same units as X
      % and Y.  PixelSizeZUnit is the X/Y pixel size in units of um.
      if isfield(DriftParams, 'PixelSizeZUnit')
         PixelSizeZUnit = DriftParams.PixelSizeZUnit;
      else
         PixelSizeZUnit = 0.1;   % default value in um
      end
      % L_intra is the intra-dataset threshold used in the cost function
      % computation to limit spurious nearest neighbor distances.  L_intra is
      % in units of pixels.
      if isfield(DriftParams, 'L_intra')
         L_intra = DriftParams.L_intra;
      else
         L_intra = 1;   % default value in pixels
      end
      % L_inter is the inter-dataset threshold used in the cost function
      % computation to limit spurious nearest neighbor distances.  L_inter is
      % in units of pixels.
      if isfield(DriftParams, 'L_inter')
         L_inter = DriftParams.L_inter;
      else
         L_inter = 2;   % default value in pixels
      end
      PDegree      = DriftParams.PDegree;
      TolFun_intra = DriftParams.TolFun_intra;
      TolX_intra   = DriftParams.TolX_intra;
      TolFun_inter = DriftParams.TolFun_inter;
      TolX_inter   = DriftParams.TolX_inter;
      Init_inter   = DriftParams.Init_inter;
      Verbose      = DriftParams.Verbose;

      if any(isfield(DriftParams, {'NDatasets', 'NFrames'}))
         SMD = ReorganizeDatasets(SMD, DriftParams);
         % SMD.NFrames has changed, so if Init_inter was set to the old value,
         % now reset it to the new value.
         if Init_inter == NFrames_C
            Init_inter = SMD.NFrames;
         end
      end
%  else
%     % Default values.
%     L_intra       = 1;     % intra-dataset threshold
%     L_inter       = 2;     % inter-dataset threshold
%     PixelSizeZUnit = 0.1;  % pixel size in um
%     PDegree       = 1;     % degree of the intra-dataset fitting polynomial
%                            % for drift correction
%     TolFun_intra = 1e-2;   % termination tolerance on the function value
%     TolX_intra   = 1e-4;   % termination tolerance on the fitting polynomial
%     TolFun_inter = 1e-2;   % termination tolerance on the function value
%     TolX_inter   = 1e-4;   % termination tolerance on the fitting polynomial
%     Init_inter   = SMD.NFrames; % initialization with respect to the previous
%                            % dataset for inter-dataset drift correction
%  end

   Statistics.Ndims          = Ndims;
   Statistics.PixelSizeZUnit = PixelSizeZUnit;
   Statistics.NDatasets      = SMD.NDatasets;
   Statistics.NFrames        = SMD.NFrames;
   Statistics.NDatasets_C    = NDatasets_C;
   Statistics.NFrames_C      = NFrames_C;
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

   SMRS = cell(1, SMD.NDatasets);
   % Count the number of iterations and function calls.
   it = 0;   fc = 0;
   % Note that variables like X, Y, Z are vectors, while variables like XY are
   % n x Ndims matrices.  The basic reason for combining vectors into matrices
   % is to allow easier generalization to 3D.
   for i = 1:SMD.NDatasets
      mask = SMD.DatasetNum == i;
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
         if exitflag ~= 1 && Verbose >= 1
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
      SMRS{i}.XY = XYC;
      SMRS{i}.n  = n;

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

   % ---------- Inter-Dataset drift correction --------------------------------

   tic;
   % These options will be used by fminsearch below.
   options = optimset('TolFun', TolFun_inter, 'TolX', TolX_inter);
   %options = optimset('TolFun', TolFun_inter, 'TolX', TolX_inter, ...
   %                  'PlotFcns', @optimplotfval);

   % Create a NeighborSearcher object for k-nearest neighbors search from the
   % first dataset.
   NS = createns(SMRS{1}.XY);
   % Count the number of iterations and function calls.
   it = 0;   fc = 0;
   for i = 2:SMD.NDatasets
      XY2 = SMRS{i}.XY;

      % Initialize the inter-dataset constant shift from the previous dataset's
      % initial or final drift depending on the value of Init_inter (or 0 for
      % no expected drift).  Note that we assume here that the inter-dataset
      % drift correction is constant.
      if Init_inter == 0
         P0_inter = zeros(1, Ndims);
      else
         if Ndims == 2
            P0_inter = [ SMD.DriftX(Init_inter, i - 1), ...
                         SMD.DriftY(Init_inter, i - 1) ];
         else
            P0_inter = [ SMD.DriftX(Init_inter, i - 1), ...
                         SMD.DriftY(Init_inter, i - 1), ...
                         SMD.DriftZ(Init_inter, i - 1) ];
         end
      end
      P0 = double(P0_inter);
      if SMRS{i}.n <= 1
         XY2C = XY2;
         P = P0;
      else
         [P, ~, exitflag, output] = ...
            fminsearch(@minD_inter, P0, options, NS, XY2, Ndims, L_inter);
         if exitflag ~= 1 && Verbose >= 1
            fprintf( ...
               'driftCorrectKNN fminsearch on minD_inter exitflag = %d\n', ...
                    exitflag);
         end
         it = it + output.iterations;
         fc = fc + output.funcCount;

         [~, XY2C] = minD_inter(P, NS, XY2, Ndims, L_inter);
      end

      % Values corrected for drift.
      SMRS{i}.XY = XY2C;

      % Update the drift correction for dataset i.
      SMD.DriftX(:, i) = SMD.DriftX(:, i) + P(1);
      SMD.DriftY(:, i) = SMD.DriftY(:, i) + P(2);

      if Ndims == 3
         SMD.DriftZ(:, i) = SMD.DriftZ(:, i) + P(3);
      end
   end
   Statistics.Inter_iterations  = it;
   Statistics.Inter_funcCount   = fc;
   Statistics.Inter_elapsedTime = toc;

   % --------------------------------------------------------------------------

   % Make sign of drift consistent with what was done before.
   SMD.DriftX = - SMD.DriftX;
   SMD.DriftY = - SMD.DriftY;
   if Ndims == 3
      SMD.DriftZ = - SMD.DriftZ .* PixelSizeZUnit;
   end

   % Save the drift corrected coordinates in SMD.
   XYout = cell(1, Ndims);
   for i = 1:SMD.NDatasets
      for j = 1:Ndims
         XYout{j} = [XYout{j}; SMRS{i}.XY(:, j)];
      end
   end
   SMD.X = XYout{1}(idx);
   SMD.Y = XYout{2}(idx);
   if Ndims == 3
      SMD.Z = XYout{3}(idx) .* PixelSizeZUnit;
   end

   % If the dataset organization has been modified (because the user chose a
   % reorganization), restore the original scheme back.
   if isfield(SMD, 'Collected')
      SMD.Internal.NDatasets  = SMD.NDatasets;
      SMD.Internal.NFrames    = SMD.NFrames;
      SMD.Internal.DatasetNum = SMD.DatasetNum;
      SMD.Internal.FrameNum   = SMD.FrameNum;
      SMD.Internal.DriftX = SMD.DriftX;
      SMD.Internal.DriftY = SMD.DriftY;
      if Ndims == 3
         SMD.Internal.DriftZ = SMD.DriftZ;
      end

      SMD.NDatasets  = SMD.Collected.NDatasets;
      SMD.NFrames    = SMD.Collected.NFrames;
      SMD.DatasetNum = SMD.Collected.DatasetNum;
      SMD.FrameNum   = SMD.Collected.FrameNum;
      SMD = rmfield(SMD, 'Collected');

      DriftX = SMD.DriftX;
      SMD.DriftX = zeros(SMD.NFrames, SMD.NDatasets, 'single');
      % The notation below reshapes the internal matrix produced for DriftX
      % into the shape expected for the original NDatasets and NFrames.
      SMD.DriftX(:) = DriftX(:);
      DriftY = SMD.DriftY;
      SMD.DriftY = zeros(SMD.NFrames, SMD.NDatasets, 'single');
      SMD.DriftY(:) = DriftY(:);
      if Ndims == 3
         DriftZ = SMD.DriftZ;
         SMD.DriftZ = zeros(SMD.NFrames, SMD.NDatasets, 'single');
         SMD.DriftZ(:) = DriftZ(:);
      end
   end

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

% -----------------------------------------------------------------------------

function [sumNND, X2] = minD_inter(Theta, NS, X2, Ndims, L_inter)
% Sum of nearest neighbor distances for inter-dataset drift correction.

   % Theta(1) = DriftX, Theta(2) = DriftY, Theta(3) = DriftZ
   for j = 1:Ndims
      X2(:, j) = X2(:, j) + Theta(j);
   end

   %[~, D] = knnsearch(NS, X2, 'K', 1);
   [~, D] = knnsearch(NS, X2);
   sumNND = sum(min(D, L_inter));

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
      sumNND = sum(min(D, L_inter)) * (1 + n/numel(D));
   end
   
end

% =============================================================================

function SMD = ReorganizeDatasets(SMD, DriftParams)
% Reorganize the collected frames into user specified dataset divisions.

   % NDatasets is the number of datasets.
   % NFrames   is the number of frames per dataset.
   % NFrames_total is the total number of frames over all datasets.
   NFrames_total = SMD.NDatasets * SMD.NFrames;
   if isfield(DriftParams, 'NDatasets') && ~isempty(DriftParams.NDatasets)
      NDatasets = DriftParams.NDatasets;
      NFrames   = NFrames_total / NDatasets;
      if mod(NFrames, 1) ~= 0
         error(['DriftParams.NDatasets (%d) does not divide evenly into\n', ...
                'the total number of frames (%d)!'], NDatasets, NFrames_total);
      end
      if isfield(DriftParams, 'NFrames') && ~isempty(DriftParams.NFrames)
         if DriftParams.NFrames ~= NFrames
         error(['DriftParams.NDatasets * DriftParams.NFrames (%d * %d) !=\n',...
                'the total number of frames (%d)!'], ...
               DriftParams.NDatasets, DriftParams.NFrames, NFrames_total);
         end
      end
   elseif isfield(DriftParams, 'NFrames') && ~isempty(DriftParams.NFrames)
      NFrames   = DriftParams.NFrames;
      NDatasets = NFrames_total / NFrames;
      if mod(NDatasets, 1) ~= 0
         error(['DriftParams.NFrames (%d) does not divide evenly into\n', ...
                'the total number of frames (%d)!'], NFrames, NFrames_total);
      end
   end

   % Compute absolute frame number as if there was only one dataset.
   FrameNumAbs = ...
      (SMD.DatasetNum - 1)*double(SMD.NFrames) + double(SMD.FrameNum);
   % Reorganize the absolute frame numbers into new dataset divisions.
   DatasetNum = ones(size(FrameNumAbs));
   FrameNum   = zeros(size(FrameNumAbs));
   for i = 1 : numel(DatasetNum)
      DatasetNum(i) = ...
         (FrameNumAbs(i) - 1 - mod(FrameNumAbs(i) - 1, NFrames)) / NFrames + 1;
      FrameNum(i)   = mod(FrameNumAbs(i) - 1, NFrames) + 1;
   end

   % Save collected dataset/frame numbering.
   SMD.Collected.NDatasets  = SMD.NDatasets;
   SMD.Collected.NFrames    = SMD.NFrames;
   SMD.Collected.DatasetNum = SMD.DatasetNum;
   SMD.Collected.FrameNum   = SMD.FrameNum;
   % Replace the old numbering scheme with the new one.
   SMD.NDatasets  = NDatasets;
   SMD.NFrames    = NFrames;
   SMD.DatasetNum = DatasetNum;
   SMD.FrameNum   = FrameNum;

end
