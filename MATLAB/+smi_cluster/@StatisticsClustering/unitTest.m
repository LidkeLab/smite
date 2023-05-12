% Test calling StatisticsClustering routines.
function success = unitTest()

success = 0;

SaveDir = smi_helpers.mkSMITETmpDir('unitTest', 'StatisticsClustering');

% Examples of how to call StatisticsClustering routines.

% --- 2D ---

SMF = smi_core.SingleMoleculeFitting();
SMF.Data.ResultsDir = SaveDir;
SC = smi_cluster.StatisticsClustering(SMF);
SC.BaseName = 'twoD';
SC.Verbose = 2;
ROI = [0, 1000, 0, 1000];   % nm

if ~exist(SMF.Data.ResultsDir, 'dir')
   mkdir(SMF.Data.ResultsDir);
end

SIM = smi_sim.SimSMLM();
kmer = 6;                % Hextets
radius_kTet = 25;        % Radius of the k-tets (nm)
PixSize = 100;           % nm in a pixel
receptorDensity = 1;     % Receptors/um^2
SIM.Rho = receptorDensity / (1000 / PixSize)^2;
SIM.simkTets(kmer, radius_kTet);

XY = 1000 * rand(2000, 2);
clear RoI

% Make RoI structure.
RoI{1}.X{1} = XY(1:2:end, 1);
RoI{1}.Y{1} = XY(2:2:end, 2);
RoI{1}.X{2} = SIM.SMD_Model.X;
RoI{1}.Y{2} = SIM.SMD_Model.Y;
RoI{1}.ROI  = ROI;

% Make SMD structures.
SMD1.X = RoI{1}.X{1};
SMD1.Y = RoI{1}.Y{1};
SMD2 = SIM.SMD_Model;
SMD2.X = SMD2.X - min(SMD2.X);
SMD2.Y = SMD2.Y - min(SMD2.Y);
SC.ROI = ROI;

fprintf('2D examples:\n\n');

% Analyze a single label dataset.
particle_types = {'A'};

fprintf('pairwiseDist ...\n');
SC.pairwiseDist(particle_types, SMD1);

% hopkins can be called this way or as further below.
%fprintf('hopkins ...\n');
%SC.hopkins(particle_types, SMD1);

fprintf('ripley ...\n');
SC.ripley(particle_types, SMD1);

% Analyze two datasets containing different labels.
particle_types = {'A', 'B'};

% Note that this makes a separate calculation for each label provided in RoI.
fprintf('hopkins ...\n');
SC.hopkins(particle_types, RoI);

% pairwiseMutualDist and bivariateRipley operate on pairs of labels.
fprintf('pairwiseMutualDist ...\n');
SC.pairwiseMutualDist(particle_types, SMD1, SMD2);

% Demonstrate alternative ways to call these routines.  This methodology
% applies to all StatisticsClustering main routines.  Note that
% pairwiseMutualDist and bivariateRipley are expecting either a single
% argument RoI structure (see smi_helpers.ROITools) or two SMD structures or
% two Nx2 (or Nx3) arrays, while the other functions are expecting an RoI
% structure or one SMD structure or one array.
fprintf('bivariateRipley ...\n');
SC.bivariateRipley(particle_types, SMD1, SMD2);
SC.BaseName = 'twoD_array';
SC.bivariateRipley(particle_types, [SMD1.X, SMD1.Y], [SMD2.X, SMD2.Y]);
SC.BaseName = 'twoD_RoI';
SC.ROI = [];
SC.bivariateRipley(particle_types, RoI);

SC.BaseName = 'twoD';
SC.ROI = ROI;

% ROIcombined series.  Here, we only have one ROI (RoI{1}).
fprintf('hopkins_ROIcombined ...\n');
SC.hopkins_ROIcombined(1, RoI);

fprintf('ripley_ROIcombined ...\n');
SC.ripley_ROIcombined(1, RoI);

fprintf('bivariateRipley_ROIcombined ...\n');
SC.bivariateRipley_ROIcombined(particle_types, 1, RoI);

fprintf('Done 2D.\n');

% --- 3D ---

SMF = smi_core.SingleMoleculeFitting();
SMF.Data.ResultsDir = SaveDir;
SC = smi_cluster.StatisticsClustering(SMF);
SC.BaseName = 'threeD';
ROI = [0, 1000, 0, 1000, 0, 1000];   % nm

if ~exist(SMF.Data.ResultsDir, 'dir')
   mkdir(SMF.Data.ResultsDir);
end

XYZ = 1000 * rand(2000, 3);

% Make RoI structure.
clear RoI
RoI{1}.X{1} = XYZ(1:2:end, 1);
RoI{1}.Y{1} = XYZ(1:2:end, 2);
RoI{1}.Z{1} = XYZ(1:2:end, 3);
RoI{1}.X{2} = XYZ(2:2:end, 1);
RoI{1}.Y{2} = XYZ(2:2:end, 2);
RoI{1}.Z{2} = XYZ(2:2:end, 3);
RoI{1}.ROI = ROI;

% Make SMD structures.
SMD1.X = RoI{1}.X{1};
SMD1.Y = RoI{1}.Y{1};
SMD1.Z = RoI{1}.Z{1};
SMD2.X = RoI{1}.X{2};
SMD2.Y = RoI{1}.Y{2};
SMD2.Z = RoI{1}.Z{2};
SC.ROI = ROI;

fprintf('3D examples:\n\n');

particle_types = {'A'};

fprintf('PairwiseDist ...\n');
SC.pairwiseDist(particle_types, SMD1);

fprintf('Hopkins ...\n');
SC.hopkins(particle_types, SMD1);

fprintf('Ripley ...\n');
SC.ripley(particle_types, SMD1);

particle_types = {'A', 'B'};

fprintf('PairwiseMutualDist ...\n');
SC.pairwiseMutualDist(particle_types, SMD1, SMD2);

fprintf('BivariateRipley ...\n');
SC.bivariateRipley(particle_types, SMD1, SMD2);

fprintf('Done 3D.\n');

fprintf('All results in %s\n', SMF.Data.ResultsDir);

success = 1;

end
