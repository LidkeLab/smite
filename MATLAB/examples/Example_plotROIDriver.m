close all

%plotROIDriver(PixelSize, options, start_datadir, SaveDir)

%options = {'MAPN', 'Gaussian', 'Boundary', 'Cluster', 'ROIImages'};
options = {'MAPN', 'Gaussian', 'Boundary', 'Cluster', 'OneImage', 'NoSave'};

start_DataDir = '/mnt/nas/cellpath/Personal Folders/Rachel/20240124_IgE-AntigenStrain';

SaveDir = '/mnt/nas/lidkelab/Personal Folders/MJW/ROI/NEW/OUT';

CI = smi_cluster.ClusterInterface;
CI.plotROIDriver(97.8, options, start_DataDir, SaveDir);
