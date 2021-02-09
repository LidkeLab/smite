classdef DriftCorrection < handle

% driftCorrectKNN performs drift correction on 2D or 3D data
% provided in an SMD structure using K nearest neighbor (KNN) searching,
% returning an updated structure with drift corrected coordinates.  Plots of
% the drift estimates can be produced with plotDriftCorrection and some
% additional measures with calcDCResidual.
%
% EXAMPLE USAGE (see also unitTest):
%    DC = smi_core.DriftCorrection(SMF, SMDin);
%    SMDIntra = [];
%    for i = 1 : NDatasets
%       [SMDIntra_i, StatisticsIntra] = DC.driftCorrectKNNIntra(SMDin_i, i);
%       SMDIntra = smi_core.SingleMoleculeData.catSMD(SMDIntra, SMDIntra_i);
%    end
%    [SMDInter, StatisticsInter] = DC.driftCorrectKNNInter(SMDIntra);
%    SMDout = SMDInter;

% =============================================================================
properties

   % Intra-dataset threshold (pixel)
   L_intra        = 1;
   % Inter-dataset threshold (pixel)
   L_inter        = 2;
   % X/Y pixel size in um (only needed for 3D drift correction)
   PixelSizeZUnit = 0.1;
   % Degree of the intra-dataset fitting polynomial for drift rate
   PDegree        = 1;
   % Termination tolerance on the intra-dataset function value
   TolFun_intra   = 1e-2;
   % Termination tolerance on the intra-dataset fitting polynomial
   TolX_intra     = 1e-4;
   % Termination tolerance on the inter-dataset function value
   TolFun_inter   = 1e-2;
   % Termination tolerance on the inter-dataset fitting polynomial
   TolX_inter     = 1e-4;
   % Initialization wrt the previous dataset for inter-dataset drift correction
   % The value should be either 0 (no initial drift), 1 (initial drift of the
   % previous dataset) or SMD.NFrames (final drift); zero or initial drift
   % should work well with brightfield registration, while final drift works
   % well generally (but the optimization process may not converge quite as
   % quickly).
   Init_inter     = 0;
   % Semi-redundant variable, needed because SMD may not exist when the
   % constructor is invoked, so Init_inter has to be set in
   % driftCorrectKNNInter (only needed when breaking intra-dataset and
   % inter-dataset calculations up).
   BFRegistration = true;
   % If non-empty, override the collected value of number of datasets
   NDatasets      = [];
   % If non-empty, override the collected value of number of frames per dataset
   NFrames        = [];
   % Verbosity level
   Verbose        = 1;

end % properties
% =============================================================================

% =============================================================================
properties(SetAccess = protected)

   % Indexing array to record rearrangements of points into datasets
   idx;
   % Values corrected for drift; fields: XY, n
   SMRS = {};

end % properties(SetAccess = protected)
% =============================================================================

% =============================================================================
methods

   [SMD, Statistics] = driftCorrectKNN(obj, SMD)
   [SMD, Statistics] = driftCorrectKNNIntra(obj, SMD, iDataset)
   [SMD, Statistics] = driftCorrectKNNInter(obj, SMD)
   DC_fig = plotDriftCorrection(obj, SMD, option)

   % Constructor.
   function obj = DriftCorrection(SMF, SMD)
   % SMF values, if provided, can override some of the class properties.
   % SMD is needed for SMD.NFrames when SMF.DriftCorrection.BFRegistration is
   % false.

      if exist('SMF', 'var')
         obj.L_intra        = SMF.DriftCorrection.L_intra;
         obj.L_inter        = SMF.DriftCorrection.L_inter;
         obj.PixelSizeZUnit = SMF.DriftCorrection.PixelSizeZUnit;
         obj.PDegree        = SMF.DriftCorrection.PDegree;
         obj.BFRegistration = SMF.DriftCorrection.BFRegistration;
         if SMF.DriftCorrection.BFRegistration
            obj.Init_inter  = 0;
         else
            if exist('SMD', 'var')
               obj.Init_inter  = SMD.NFrames;
            else
               %error('SMD not available when BFRegistration is false.');
            end
         end
      end

   end

end % methods
% =============================================================================

% =============================================================================
methods(Static)

   [residual, dist, rmse, nnfig] = calcDCResidual(SMD, X_True, Y_True, Z_True)
   [success, SMD2, SMD3, Statistics2, Statistics3] = unitTest()

end % methods(Static)
% =============================================================================

end % classdef DriftCorrection
