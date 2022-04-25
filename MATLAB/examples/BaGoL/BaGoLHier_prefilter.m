% Example pf BaGoL pre-filtering, starting with a Picasso .h5 file or a
% Results.mat file.  The input SR data should be run with no thresholding on
% max Photons, no frame connection and no drift correction (drift correction
% may be OK in some circumstances, but has not been tested).  The flow of data
% for this analysis is as follows:
%
%    SR data -> intensity filter
%            -> apply SE_Adjust
%            -> frame connection
%            -> remove connections which involve only one frame
%            -> producing a dataset like:
%               EGFR10min_IFilter_FC_SEAdjust0p07pixels.mat
%            -> NND filter (via BaGoLHier_wrapper/BaGoLHier_analysis)
%            -> BaGoL
%
% Intensity filter defaultly removes localizations that have photon counts
%    > 2 * mean(SMD.Photons)
% SE_Adjust inflates the localization standard errors by a small amount to make
%    it easier for BaGoL to group them together.
% Frame connection is standard and a dataset is produced in which each
%    localization is represented by two or more frames.
% Afterwards, BaGoL is run (via the scripts BaGoLHier_wrapper/analysis), first
%    performing a nearest neighbor distance (NND) filtering in which
%    localizations are removed if their N_NN (minimum number of nearest
%    neighbors) distances exceed 3 standard deviations from the median of the
%    localization NNDs.  N_NN is set in BaGoLHier_wrapper.  Note that this step
%    can sometimes run out of memory if the number of localizations is large.
%    Setting N_NN to zero will skip this thresholding.
%
%    Once the above is done, regular BaGoL is invoked, producing a variety of
%    files as detailed in the other scripts.

DataDir = '.';
BaseName = 'Data_2022-4-12-19-51-9_SRZoom20_Results';

SaveDir = DataDir;
FileName = [BaseName, '.mat'];
%SMD = smi.BaGoL.loadPICASSOh5(DataDir, FileName);
load(fullfile(DataDir, FileName));
SMDCopy = SMD;

%% Filter localizations based on intensity.
MeanMultiplier = 2; % intensity > MeanMultiplier*mean(Intensity) are removed
MaxPhotons = MeanMultiplier * mean(SMD.Photons);
SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, SMD.Photons<=MaxPhotons);

%% Inflate standard errors.

SEAdjust = 0.07; % pixels, make 0.\d\d so saving below displays correctly
SMD.X_SE = SMD.X_SE + SEAdjust;
SMD.Y_SE = SMD.Y_SE + SEAdjust;

%% Perform frame connection
N = numel(SMD.X);
%SMD.DatasetNum = ones(size(SMD.FrameNum), 'uint32');
%SMD.ThreshFlag = zeros(size(SMD.FrameNum));
%SMD.Photons = zeros(size(SMD.FrameNum));
%SMD.Bg = zeros(size(SMD.FrameNum));
%SMD.LogLikelihood = zeros(size(SMD.FrameNum));
%SMD.NDims = 2;
%SMD.NDatasets = 10;
%SMD.NFrames = max(SMD.FrameNum);
%SMD.XSize = 2 ^ nextpow2(max([SMD.X; SMD.Y]));
%SMD.YSize = SMD.XSize;
SMF = smi_core.SingleMoleculeFitting;
SMF.FrameConnection.Method = 'Hypothesis test';
SMF.FrameConnection.LoS = 0.01;
SMF.FrameConnection.MaxSeparation = 1; % pixels
SMF.FrameConnection.MaxFrameGap = 5; % frames
FC = smi_core.FrameConnection(SMD, SMF);
FC.performFrameConnection();
SMD = FC.SMDCombined;

%% Remove localizations representing fewer than 2 frame-connected localizations.
SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, SMD.NCombined > 1);
save(fullfile(DataDir, sprintf('%s_IFilter_FC_SEAdjust0p%02ipixels.mat', ...
    BaseName, round(SEAdjust*100))), ...
    'SMD');
