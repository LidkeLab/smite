classdef DriftCorrection < handle

% driftCorrectKNN performs drift correction on 2D or 3D data
% provided in an SMD structure using K nearest neighbor (KNN) searching,
% returning an updated structure with drift corrected coordinates.  Plots of
% the drift estimates can be produced with plotDriftCorrection and some
% additional measures with calcDCResidual.

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

end % properties
% =============================================================================

% =============================================================================
methods

   [SMD, Statistics] = driftCorrectKNN(SMD, DriftParams)
   DC_fig = plotDriftCorrection(SMD, DriftParams, option)

   % Constructor.
   function obj = DriftCorrection(obj, SMF)
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
   [SMD2, SMD3, Statistics2, Statistics3] = unitTest()

end % methods(Static)
% =============================================================================

end % classdef DriftCorrection
