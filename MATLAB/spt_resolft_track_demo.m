

addpath('mex\')
addpath('ptx\')
%% here only need to select the data, other parameters will be set below
SMF = smi_core.SingleMoleculeFitting;
SMF.gui

%% set path to dark frame file
offsetfile = 'Y:\sCMOS Calibrations\SPT\GainCalibration_darkFrames_2022_10_26_18_02_10.mat';
load(offsetfile)
ccdoffset = mean(sequence,3);
ccdvar = var(single(sequence),1,3);

%% Prepare the SPT class object.

SPT = smi.SPT(SMF, false);
SPT.GenerateMovies = false;
SPT.GeneratePlots = false;
SPT.SMF.Data.CameraOffset = ccdoffset;
SPT.SMF.Data.CameraGain = 3.8*ones(size(ccdoffset));
SPT.SMF.Data.CameraReadNoise = ccdvar;
SPT.SMF.Data.PixelSize = 0.1; % um
SPT.SMF.BoxFinding.MinPhotons = 20;
SPT.SMF.BoxFinding.BoxSize = 11;
SPT.SMF.Fitting.FitType = 'XYNBSXSY';
SPT.SMF.Fitting.NParams = 6;
SPT.SMF.Fitting.PSFSigma = 2;
SPT.SMF.Fitting.Iterations = 50;
SPT.SMF.Thresholding.AutoThreshLogL = 1;
SPT.SMF.Thresholding.MaxXY_SE = 1;
SPT.SMF.Thresholding.MaxPSFSigma = 3;
SPT.SMF.Thresholding.MinPhotons = 10;
SPT.SMF.Thresholding.MinPValue = 0;
SPT.SMF.FrameConnection.On = 0;
SPT.SMF.DriftCorrection.On = 0;
SPT.performFullAnalysis();

%% show tracks > 20 localizations
Nf = length(SPT.TR);
figure;
ha = axes;hold on
for tt=1:Nf
    trk = SPT.TR(tt);
    if numel(trk.X)>20

        plot3(trk.X,trk.Y,trk.FrameNum,'.-')
    end
end
xlabel('x (pixel)')
ylabel('y (pixel)')
zlabel('frame')
view(30,30)
grid(ha,"on")
%% run change point analysis for tracks > 20 localizations
h = figure('Position',[200,100,600,500]);
tg1 = uitabgroup(h); % tabgroup
LogBayesThreshold = 150;
trackID = []; % save tracks with change points
count = 1;
for ii = 1:Nf
    trk = SPT.TR(ii);
    if numel(trk.X)>20
        icp = smi_stat.ChangeDetection(round(double(trk.Photons)),LogBayesThreshold);
        if ~isempty(icp.ChangePoints)
            xs=1:icp.Nobservations;
            tab1 = uitab(tg1,'title',num2str(ii));
            ha = axes('parent',tab1);
            stairs(xs, icp.Data, '-ko');
            hold on;
            Is=smi_stat.ChangeDetection.modelIntensity(icp.Nobservations, icp.ChangePoints, icp.Intensity);
            stairs(xs, Is, '-r', 'LineWidth', 2.0);
            yl=ylim;
            ylim([0,yl(2)]);
            hold off;
            xlabel('time');
            ylabel('intensity');
            legend('Data', 'Estimated', 'Location', 'North');
            trackID(count) = ii;
            count = count+1;
        end
    end
end

%% add final bleaching step for selected candidate track from above plot
TR_candidate = [];
data = SPT.ScaledData;
h = figure('Position',[200,100,600,500]);
tg1 = uitabgroup(h); % tabgroup
for ii = 1:numel(trackID)
    ind = trackID(ii);
    trk = SPT.TR(ind);
    TR_candidate{ii} = trk;
    % crop data from last bleaching step, 10 frame
    bxsz = 11;
    roi = data(round(trk.Y(end)-bxsz/2):round(trk.Y(end)+bxsz/2),round(trk.X(end)-bxsz/2):round(trk.X(end)+bxsz/2),trk.FrameNum(end)+1:min([trk.FrameNum(end)+10,size(data,3)]));
    roi = movmean(roi,3,3);
    % localization
    NP=6;
    SZ = size(roi,1);
    BSZ=128;
    KernelID='_XYNBSXSY_';
    k = parallel.gpu.CUDAKernel('smi_cuda_gaussMLEv2.ptx','smi_cuda_gaussMLEv2.cu',KernelID);
    NFitsActual = size(roi,3);
    PSFSigma = 2;
    Iterations = 50;
    k.GridSize = [ceil(NFitsActual/BSZ) 1];
    k.ThreadBlockSize = [BSZ 1];
    d_Parameters=zeros(NFitsActual,NP,'single');
    d_CRLBs=zeros(NFitsActual,NP,'single');
    d_LogLikelihood=zeros(NFitsActual,1,'single');

    [P, CRLB,LL] = feval(k,single(roi),PSFSigma,SZ,Iterations,d_Parameters,d_CRLBs,d_LogLikelihood,NFitsActual);
    % add photons from bleaching step to original photon sequence
    photon = gather(P(:,3));
    mask = photon<300 & photon>2;
    photon = photon(mask);
    Photons_laststep =  round(double([trk.Photons;photon]));
    TR_candidate{ii}.Photons_laststep = Photons_laststep;
    % get change point for new photon trace
    icp = smi_stat.ChangeDetection(Photons_laststep,LogBayesThreshold);

    %show plots
    xs=1:icp.Nobservations;
    tab1 = uitab(tg1,'title',num2str(ind));
    ha = axes('parent',tab1);
    stairs(xs, icp.Data, '-ko');
    hold on;
    Is=smi_stat.ChangeDetection.modelIntensity(icp.Nobservations, icp.ChangePoints, icp.Intensity);
    stairs(xs, Is, '-r', 'LineWidth', 2.0);
    yl=ylim;
    ylim([0,yl(2)]);
    hold off;
    xlabel('time');
    ylabel('intensity');
    legend('Data', 'Estimated', 'Location', 'North');

end
%% overlay raw data with selected track
data = SPT.ScaledData;
ind = 1;
trk = SPT.TR(ind);
imt = zeros(size(data));
for ii = 1:length(trk.X)
imt(round(trk.Y(ii)),round(trk.X(ii)),trk.FrameNum(ii)) = 1;
end
overlay(data/median(max(data,[],[1,2]))*100,imt)

