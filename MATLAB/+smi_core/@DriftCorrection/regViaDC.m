%regViaDC (registration via drift correction) takes two differently labeled
% data collections of the same biological phenomenon and attempts to align them
% using inter-dataset drift correction.
function [delta12, Statistics] = regViaDC(SMD1, SMD2)
%
% INPUTS:
%    SMD1    fixed  data (X and Y fields in pixels)
%    SMD2    moving data (X and Y fields in pixels)
%            Data are from either BaGoL ResultsStruct files (SMR) or smite
%            Results (SMD + SMF) files taken under the same imaging conditions.
%
% OUTPUTS:
%    delta12      computed registration shift in the same units as SMD1/2.X/Y
%    Statistics   statistical information about the algorithm performance
%                 (see driftCorrectKNN)
%
% Note that:
%    drifted coordinates - drift correction = drift corrected coordinates

% Created by
%    Michael J. Wester, 2021

   DC = smi_core.DriftCorrection;

   % Intra-dataset threshold (pixel)
   DC.L_intra        = 1;
   % Inter-dataset threshold (pixel)
   DC.L_inter        = 2;
   % X/Y pixel size in um (only needed for 3D drift correction)
   DC.PixelSizeZUnit = 0.1;
   % Degree of the intra-dataset fitting polynomial for drift rate
   DC.PDegree        = 1;
   % Termination tolerance on the intra-dataset function value
   DC.TolFun_intra   = 1e-2;
   % Termination tolerance on the intra-dataset fitting polynomial
   DC.TolX_intra     = 1e-4;
   % Termination tolerance on the inter-dataset function value
   DC.TolFun_inter   = 1e-2;
   % Termination tolerance on the inter-dataset fitting polynomial
   DC.TolX_inter     = 1e-4;
   % Initialization wrt the previous dataset for inter-dataset drift correction
   % The value should be either 0 (no initial drift), 1 (initial drift of the
   % previous dataset) or SMD.NFrames (final drift); zero or initial drift
   % should work well with brightfield registration, while final drift works
   % well generally (but the optimization process may not converge quite as
   % quickly).
   DC.Init_inter     = 1;
   % If non-empty, override the collected value of number of datasets
   DC.NDatasets      = [];
   % If non-empty, override the collected value of number of frames per dataset
   DC.NFrames        = [];
   % Verbosity level
   DC.Verbose        = 1;

   % Treat the two collections as two datasets in one collection.
   SMDin.X = [SMD1.X; SMD2.X];
   SMDin.Y = [SMD1.Y; SMD2.Y];
   SMDin.NDatasets = 2;
   SMDin.NFrames   = 1;
   SMDin.DatasetNum = [repmat(1, numel(SMD1.X), 1);
                       repmat(2, numel(SMD2.X), 1)];
   SMDin.FrameNum   = repmat(1, numel(SMDin.X), 1);

   [SMDout, Statistics] = DC.driftCorrectKNN(SMDin);
   delta12 = [SMDout.DriftX(1, 2), SMDout.DriftY(1, 2)];

end
