function [SMDout, Statistics] = driftCorrectKNNInterPair(obj, SMD1, SMD2)
%driftCorrectKNNInterPair calculates inter-dataset drift directly from X,Y{,Z}
% coordinates (i.e., constant shifts between datasets).  Fitting is done via
% performing fminsearch on the (weighted) sums of nearest neighbor distances.
% NOTE: SMD1 is considered fixed in space, while SMD2 is moving.  SMD1 and SMD2
% should have the same number of datasets and frames per dataset.
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
%   SMD1, SMD2:    A structure with fields:
%      X              x coordinates (Nx1) where N is total number of points
%      Y              y coordinates (Nx1)
%      Z              z coordinates (Nx1) [OPTIONAL]
%                  Note: X and Y are typically in pixels, Z in um
%      DatasetNum     dataset number from which localization originates (Nx1)
%      FrameNum       frame   number from which localization originates (Nx1)
%      NDatasets      number of datasets
%      NFrames        number of frames in each dataset
%      DriftX         previously computed x drift (NFrames x NDatasets)
%      DriftY         previously computed y drift (NFrames x NDatasets)
%      DriftZ         previously computed z drift (NFrames x NDatasets)
%                     [OPTIONAL]
%   obj:           [class properties]
%                     optimization parameters with the following fields:
%      L_inter        inter-dataset threshold (Default = 2 pixels)
%      PixelSizeZUnit pixel size in um (units of Z; only needed for 3D)
%                     (Default = 0.1)
%      TolFun_inter   termination tolerance on the inter-dataset function value
%                     (Default = 1e-2)
%      TolX_inter     termination tolerance on the inter-dataset fitting
%                     polynomial (Default = 1e-4)
%      P0_inter       inter-dataset initialization
%                     (Default = [] for all zeros)
%      Verbose        verbosity level (Default = 1)
%
% OUTPUTS:
%   SMDout:      SMD data structure with updated fields:
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
%      Inter_iterations   inter-dataset number of fminsearch iterations
%      Inter_funcCount    inter-dataset number of fminsearch function evals
%      Inter_cost         inter-dataset cost at the found minimum (sumNND)
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

   if isfield(SMD1, 'Z') && numel(SMD1.Z) == numel(SMD1.X)
      Ndims = 3;
   else
      Ndims = 2;
   end

   DriftParams.L_inter        = obj.L_inter;
   DriftParams.PixelSizeZUnit = obj.PixelSizeZUnit;
   DriftParams.TolFun_inter   = obj.TolFun_inter;
   DriftParams.TolX_inter     = obj.TolX_inter;
   DriftParams.P0_inter       = obj.P0_inter;
   DriftParams.Verbose        = obj.Verbose;

   % Initialize various parameters.
   % PixelSizeZUnit is needed for 3D to convert Z into the same units as X
   % and Y.  PixelSizeZUnit is the X/Y pixel size in units of um.
   PixelSizeZUnit = DriftParams.PixelSizeZUnit;
   % L_inter is the inter-dataset threshold used in the cost function
   % computation to limit spurious nearest neighbor distances.  L_inter is
   % in units of pixels.
   L_inter = DriftParams.L_inter;
   TolFun_inter = DriftParams.TolFun_inter;
   TolX_inter   = DriftParams.TolX_inter;
   P0_inter     = DriftParams.P0_inter;
   Verbose      = DriftParams.Verbose;

   Statistics.Ndims          = Ndims;
   Statistics.PixelSizeZUnit = PixelSizeZUnit;
   Statistics.NDatasets      = SMD1.NDatasets;
   Statistics.NFrames        = SMD1.NFrames;
   Statistics.L_inter        = L_inter;
   Statistics.TolFun_inter   = TolFun_inter;
   Statistics.TolX_inter     = TolX_inter;

   % ---------- Intra-Dataset drift correction --------------------------------

   % ---------- Inter-Dataset drift correction --------------------------------

   %SMRS = obj.SMRS;
   % Sum image of SMD over all datasets in SMD.
   if Ndims == 2
      SMRS{1}.XY = [SMD1.X, SMD1.Y];
   else
      SMRS{1}.XY = [SMD1.X, SMD1.Y, SMD1.Z ./ PixelSizeZunit];
   end
   SMRS{1}.n  = size(SMRS{1}.XY, 1);
   % sum image of SMD2 over all datasets in SMD2.
   if Ndims == 2
      SMRS{2}.XY = [SMD2.X, SMD2.Y];
   else
      SMRS{2}.XY = [SMD2.X, SMD2.Y, SMD2.Z ./ PixelSizeZunit];
   end
   SMRS{2}.n  = size(SMRS{2}.XY, 1);
   DriftX = zeros(SMD1.NFrames, 1, 'single');
   DriftY = zeros(SMD1.NFrames, 1, 'single');
   if Ndims == 3
      DriftZ = zeros(SMD1.NFrames, 1);
   end

   tic;
   % These options will be used by fminsearch below.
   options = optimset('TolFun', TolFun_inter, 'TolX', TolX_inter);
   %options = optimset('TolFun', TolFun_inter, 'TolX', TolX_inter, ...
   %                  'PlotFcns', @optimplotfval);

   % Create a NeighborSearcher object for k-nearest neighbors search from the
   % first dataset.
   if SMRS{1}.n == 0
      % This whole procedure fails if the first dataset is empty!
      error(['First dataset is empty!  ', ...
             'Cannot perform inter-dataset drift correction.']);
   end
   NS = createns(SMRS{1}.XY);
   cost_inter = 0;
   % Count the number of iterations and function calls.
   it = 0;   fc = 0;
%  for i = 2:SMD1.NDatasets
      XY2 = SMRS{2}.XY;

      % Initialize the inter-dataset constant shift.  Note that we assume here
      % that the inter-dataset drift correction is constant.
      if isempty(P0_inter)
         P0_inter = zeros(1, Ndims);
      end
      P0 = double(P0_inter);
      if SMRS{2}.n <= 1
         XY2C = XY2;
         P = P0;
      else
         [P, ~, exitflag, output] = ...
            fminsearch(@smi_core.DriftCorrection.minD, ...
                       P0, options, XY2, [], Ndims, L_inter, NS);
         if exitflag ~= 1 && Verbose >= 1
            fprintf( ...
               'driftCorrectKNN fminsearch on minD_inter exitflag = %d\n', ...
                    exitflag);
         end
         it = it + output.iterations;
         fc = fc + output.funcCount;

         [cost_inter, XY2C] = ...
            smi_core.DriftCorrection.minD(P, XY2, [], Ndims, L_inter, NS);
      end

      % Values corrected for drift.
      SMRS{2}.XY = XY2C;

      % Update the drift correction.
      DriftX(:, 1) = DriftX(:, 1) + P(1);
      DriftY(:, 1) = DriftY(:, 1) + P(2);

      if Ndims == 3
%        SMD1.DriftZ(:, i) = SMD1.DriftZ(:, i) + P(3);
         DriftZ(:, 1) = DriftZ(:, 1) + P(3);
      end
%  end
   Statistics.Inter_iterations  = it;
   Statistics.Inter_funcCount   = fc;
   Statistics.Inter_cost        = cost_inter;
   Statistics.Inter_elapsedTime = toc;

   % --------------------------------------------------------------------------

   SMDout = SMD1;
   SMDout.DriftX = DriftX;
   SMDout.DriftY = DriftY;
   if Ndims == 3
      SMDout.DriftZ = DriftZ .* PixelSizeZUnit;
   end

   % No longer needed when sign in minD for X updates is reversed.
%  % Make sign of drift consistent with what was done before.
%  SMDout.DriftX = - SMDout.DriftX;
%  SMDout.DriftY = - SMDout.DriftY;
%  if Ndims == 3
%     SMDout.DriftZ = - SMDout.DriftZ .* PixelSizeZUnit;
%  end

   % Save the drift corrected coordinates in SMDout.
   XYout = cell(1, Ndims);
%  for i = 1:SMD1.NDatasets
      for j = 1:Ndims
         XYout{j} = [XYout{j}; SMRS{2}.XY(:, j)];
      end
%  end
   SMDout.X = XYout{1};
   SMDout.Y = XYout{2};
   if Ndims == 3
      SMDout.Z = XYout{3} .* PixelSizeZUnit;
   end

end
