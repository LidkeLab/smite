classdef DriftCorrection < handle

% driftCorrectKNN performs drift correction on 2D or 3D data
% provided in an SMD structure using K nearest neighbor (KNN) searching,
% returning an updated structure with drift corrected coordinates.  Plots of
% the drift estimates can be produced with plotDriftCorrection and some
% additional measures with calcDCResidual.
%
% EXAMPLE USAGE (see also unitTest):
%    DC = smi_core.DriftCorrection(SMF);
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
   % degree of the intra-dataset fitting polynomial for drift rate
   PDegree        = 1;
   % termination tolerance on the intra-dataset function value
   TolFun_intra   = 1e-2;
   % termination tolerance on the intra-dataset fitting polynomial
   TolX_intra     = 1e-4;
   % termination tolerance on the inter-dataset function value
   TolFun_inter   = 1e-2;
   % termination tolerance on the inter-dataset fitting polynomial
   TolX_inter     = 1e-4;
   % initialization wrt the previous dataset for inter-dataset drift correction
   Init_inter     = 0;
   % if non-empty, override the collected value of number of datasets
   NDatasets      = [];
   % if non-empty, override the collected value of number of frames per dataset
   NFrames        = [];

end % properties
% =============================================================================

% =============================================================================
properties(SetAccess = protected)

   % indexing array to record rearrangements of points into datasets
   idx;
   % values corrected for drift; fields: X, Y, n
   SMRS = {};

end % properties(SetAccess = protected)
% =============================================================================

% =============================================================================
methods

   [SMD, Statistics] = driftCorrectKNN(obj, SMD)
   [SMD, Statistics] = driftCorrectKNNIntra(obj, SMD, iDataset)
   [SMD, Statistics] = driftCorrectKNNInter(obj, SMD)
   SCobj = initializeDriftCorrection(SCobj)
   DC_fig = plotDriftCorrection(obj, SMD, option)

   % Constructor.
   function obj = DriftCorrection(SMF)
   % SMF values, if provided, can override some of the class properties.

      if exist('SMF', 'var')
         obj.L_intra        = SMF.DriftCorrection.L_intra;
         obj.L_inter        = SMF.DriftCorrection.L_inter;
         obj.PixelSizeZUnit = SMF.DriftCorrection.PixelSizeZUnit;
         obj.PDegree        = SMF.DriftCorrection.PDegree;
         obj.Init_inter     = SMF.DriftCorrection.Init_inter;
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
