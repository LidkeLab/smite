function [SMD, Statistics] = driftCorrectKNNInter(obj, SMD)
%driftCorrectKNNInter calculates inter-dataset drift directly from X,Y{,Z} coordinates
% by fitting a polynomial depending on time (i.e., frame number) to the frames
% with each dataset (intra-dataset), and fitting constant shifts between
% datasets (inter-dataset).  Fitting is done via performing fminsearch on the
% (weighted) sums of nearest neighbor distances.  Inter-dataset portion.
%
% Sign convention:
%
%    N = numel(SMD_drifted.X);
%    for k = 1:N
%       i = SMD_corrected.FrameNum(k);
%       j = SMD_corrected.DatasetNum(k);
%       SMD_corrected.X(k) = SMD_drifted.X(k) - SMD_corrected.DriftX(i, j);
%       SMD_corrected.Y(k) = SMD_drifted.Y(k) - SMD_corrected.DriftY(i, j);
%    end
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
%      DriftX         intra-dataset x drift (NFrames x NDatasets)
%      DriftY         intra-dataset y drift (NFrames x NDatasets)
%      DriftZ         intra-dataset z drift (NFrames x NDatasets) [OPTIONAL]
%                  Note: DriftX/Y/Z should be computed by driftCorrectKNNIntra
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
%      BFRegistration override value for Init_inter from obj if BFRegistration
%                     is false and the current value for Init_inter is zero, so
%                     the two are in conflict.  BFRegistration will take
%                     precedence
%      Verbose        verbosity level (Default = 1)
%
% OUTPUTS:
%   SMD:         SMD data structure with updated fields:
%      X              drift corrected x coordinates (Nx1)
%      Y              drift corrected y coordinates (Nx1)
%      Z              drift corrected z coordinates (Nx1) [OPTIONAL]
%      DriftX         found x drift per frame (NFrames x NDatasets)
%      DriftY         found y drift per frame (NFrames x NDatasets)
%      DriftZ         found z drift per frame (NFrames x NDatasets) [OPTIONAL]
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
%   NOTE: SMD.DriftX/Y/Z are the drift corrections defined such that
%         drifted coordinates - drift correction = drift corrected coordinates
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
   DriftParams.Verbose        = obj.Verbose;
   % Override value for Init_inter from obj if BFRegistration is false and the
   % current value for Init_inter is zero, so the two are in conflict.
   % BFRegistration will take precedence.
   if ~obj.BFRegistration && DriftParams.Init_inter == 0
      DriftParams.Init_inter = SMD.NFrames; 
   end
   if ~isempty(obj.NDatasets)
      DriftParams.NDatasets   = obj.NDatasets;
      obj.Init_inter          = obj.NDatasets;
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
   Verbose      = DriftParams.Verbose;

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

   % ---------- Intra-Dataset drift correction --------------------------------

   % ---------- Inter-Dataset drift correction --------------------------------

   idx  = obj.idx;
   SMRS = obj.SMRS;

   tic;
   % These options will be used by fminsearch below.
   options = optimset('TolFun', TolFun_inter, 'TolX', TolX_inter);
   %options = optimset('TolFun', TolFun_inter, 'TolX', TolX_inter, ...
   %                  'PlotFcns', @optimplotfval);

   % Create a NeighborSearcher object for k-nearest neighbors search from the
   % first dataset.
   if SMD.NDatasets > 1 && SMRS{1}.n == 0
      % This whole procedure fails if the first dataset is empty!
      error(['First dataset is empty!  ', ...
             'Cannot perform inter-dataset drift correction.']);
   end
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
