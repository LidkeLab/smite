% An example I used (MJW) to develop SMLM.
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir           = ...
   'Y:\Sandeep\Genmab\10082020\Wien133_LQT_CD52_HexElect1\Cell_01\Label_01';
SMF.Data.FileName         = {'Data_2020-10-8-17-58-54.h5'};
SMF.Data.ResultsDir       = 'Y:\MJW\SR\Results';
SMF.Data.CameraType       = 'EMCCD';
SMF.Data.CameraGain       = 1;
SMF.Data.CameraOffset     = 0;
SMF.Thresholding.On       = true;
SMF.Thresholding.MaxXY_SE = 0.1;
SMF.FrameConnection.On    = true;
SMF.DriftCorrection.On    = true;
SMLMobj = smi.SMLM(SMF);
SMLMobj.Verbose = 1;
SMLMobj.fullAnalysis();
