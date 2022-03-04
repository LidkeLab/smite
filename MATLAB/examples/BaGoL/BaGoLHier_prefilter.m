% Example pf BaGoL pre-filtering, starting with a Picasso .h5 file.
%
%    SR data -> intensity filter
%            -> apply SE_Adjust
%            -> frame connection
%            -> remove connections which involve only one frame
%            -> producing a dataset like:
%               EGFR10min_IFilter_FC_SEAdjust0p05pixels.mat
%            -> NND filter (via BaGoL_wrapper/BaGoL_analysis)
%            -> BaGoL
%
% DataDir = '\\rayleigh.phys.unm.edu\data\BaGoL\Jungmann\kindlin_data';
DataDir = '/mnt/nas/lidkelab/BaGoL/Jungmann/kindlin_data';
SaveDir = DataDir;
FileName = 'ta_ki_int_R1_200pM_561_130mW_v1_1_MMStack_Pos0.ome_locs_render_RCC500.hdf5';
SMD = smi.BaGoL.loadPICASSOh5(DataDir, FileName);
BaseName = 'ta_ki_int_R1_200pM';
save(fullfile(DataDir, [BaseName, '_NoFilterNoFC.mat']), 'SMD')
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
SMD.FrameNum = SMD.FrameNum + 1; % starts at 0 in the above file!
SMD.DatasetNum = ones(size(SMD.FrameNum), 'uint32');
SMD.ThreshFlag = zeros(size(SMD.FrameNum));
SMD.Photons = zeros(size(SMD.FrameNum));
SMD.Bg = zeros(size(SMD.FrameNum));
SMD.LogLikelihood = zeros(size(SMD.FrameNum));
SMD.NDims = 2;
SMD.NDatasets = 1;
SMD.NFrames = max(SMD.FrameNum);
SMD.XSize = 2 ^ nextpow2(max([SMD.X; SMD.Y]));
SMD.YSize = SMD.XSize;
SMF = smi_core.SingleMoleculeFitting;
SMF.FrameConnection.Method = 'Hypothesis test';
SMF.FrameConnection.LoS = 0.01;
SMF.FrameConnection.MaxSeparation = 1; % pixels
SMF.FrameConnection.MaxFrameGap = 5; % frames
FC = smi_core.FrameConnection(SMD, SMF);
FC.performFrameConnection()
SMD = FC.SMDCombined;

%% Remove localizations representing fewer than 2 frame-connected localizations.
SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, SMD.NCombined > 1);
save(fullfile(DataDir, sprintf('%s_IFilter_FC_SEAdjust0p%02ipixels.mat', ...
    BaseName, round(SEAdjust*100))), ...
    'SMD')
