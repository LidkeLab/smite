### smite namespaces and classes

|***SMITE classes***||
-------------|---
[MATLAB/+smi](../MATLAB/+smi/README.md):
&nbsp;&nbsp;&nbsp;[@BaGoL/](../MATLAB/+smi/@BaGoL/README.md)                    | Bayesian Grouping of Localizations
&nbsp;&nbsp;&nbsp;[@Publish/](../MATLAB/+smi/@Publish/README.md)                | batch processing of SR data
&nbsp;&nbsp;&nbsp;[@SMLM/](../MATLAB/+smi/@SMLM/README.md)                      | single molecule localization microscopy
&nbsp;&nbsp;&nbsp;[@SPT/](../MATLAB/+smi/@SPT/README.md)                        | single-particle tracking analysis
[MATLAB/+smi_cluster](../MATLAB/+smi_cluster/README.md):
&nbsp;&nbsp;&nbsp;[@Clustering/](../MATLAB/+smi_cluster/@Clustering/README.md)  | clustering algorithms (DBSCAN, Voronoi, H-SET)
&nbsp;&nbsp;&nbsp;[@ClusterInterface/](../MATLAB/+smi_cluster/@ClusterInterface/README.md)         | interface functions for Clustering class
&nbsp;&nbsp;&nbsp;[@PairAnalysis/](../MATLAB/+smi_cluster/@PairAnalysis/README.md)                 | interface functions for PairCorrelation class
&nbsp;&nbsp;&nbsp;[@PairCorrelation/](../MATLAB/+smi_cluster/@PairCorrelation/README.md)           | pair correlation (auto- and cross-) on ROIs
&nbsp;&nbsp;&nbsp;[@StatisticsClustering/](../MATLAB/+smi_cluster/@StatisticsClustering/README.md) | clustering statistical analyses
&nbsp;&nbsp;&nbsp;[MATLAB/+smi_core](../MATLAB/+smi_core/README.md):
&nbsp;&nbsp;&nbsp;[@ChannelRegistration/](../MATLAB/+smi_core/@ChannelRegistration/README.md)      | channel registration
&nbsp;&nbsp;&nbsp;[@DataToPhotons/](../MATLAB/+smi_core/@DataToPhotons/README.md)                  | convert raw data to photons
&nbsp;&nbsp;&nbsp;[@DriftCorrection/](../MATLAB/+smi_core/@DriftCorrection/README.md)              | drift correction on 2D/3D data
&nbsp;&nbsp;&nbsp;[@FRC/](../MATLAB/+smi_core/@FRC/README.md)                   | Fourier Ring Correlation for image resolution
&nbsp;&nbsp;&nbsp;[@FindROI/](../MATLAB/+smi_core/@FindROI/README.md)           | find/collate subregions from a 2D image stack
&nbsp;&nbsp;&nbsp;[@FrameConnection/](../MATLAB/+smi_core/@FrameConnection/README.md)              | frame connection
&nbsp;&nbsp;&nbsp;[@LoadData/](../MATLAB/+smi_core/@LoadData/README.md)         | load raw microscope data from .mat/.h5 files
&nbsp;&nbsp;&nbsp;[@LocalizeData/](../MATLAB/+smi_core/@LocalizeData/README.md) | find localizations in raw data
&nbsp;&nbsp;&nbsp;[@SingleMoleculeData/](../MATLAB/+smi_core/@SingleMoleculeData/README.md)        | define SMD structure
&nbsp;&nbsp;&nbsp;[@SingleMoleculeFitting/](../MATLAB/+smi_core/@SingleMoleculeFitting/README.md)  | define SMF structure
&nbsp;&nbsp;&nbsp;[@Threshold/](../MATLAB/+smi_core/@Threshold/README.md)       | threshold based on various SMLM properties
&nbsp;&nbsp;&nbsp;[@TrackingResults/](../MATLAB/+smi_core/@TrackingResults/README.md)              | define Tracking Results (TR) structure
&nbsp;&nbsp;&nbsp;[GaussMLE.m](../MATLAB/+smi_core/GaussMLE.m/README.md)        | max likelihood estimate of 2D Gaussian blobs
[MATLAB/+smi_helpers](../MATLAB/+smi_helpers):
&nbsp;&nbsp;&nbsp;[@Filters/](../MATLAB/+smi_helpers/@Filters/README.md)        | filters useful for BaGoL operating on SMDs
&nbsp;&nbsp;&nbsp;[@ROITools/](../MATLAB/+smi_helpers/@ROITools/README.md)      | select ROIs from an image; save in a structure
[MATLAB/+smi_psf](../MATLAB/+smi_psf):
&nbsp;&nbsp;&nbsp;[@PointSpreadFunction/](../MATLAB/+smi_psf/@PointSpreadFunction/README.md)       | create and quantify point spread functions
&nbsp;&nbsp;&nbsp;[@Zernike/](../MATLAB/+smi_psf/@Zernike/README.md)            | low-level Zernike polynomial functions
[MATLAB/+smi_sim](../MATLAB/+smi_sim):
&nbsp;&nbsp;&nbsp;[@GaussBlobs/](../MATLAB/+smi_sim/@GaussBlobs/README.md)      | generate 2D Gaussian blob images
&nbsp;&nbsp;&nbsp;[@SimSMLM/](../MATLAB/+smi_sim/@SimSMLM/README.md)            | simulate SMLM data
&nbsp;&nbsp;&nbsp;[@SimSPT/](../MATLAB/+smi_sim/@SimSPT/README.md)              | simulate SPT data
[MATLAB/+smi_stat](../MATLAB/+smi_stat/README.md):
&nbsp;&nbsp;&nbsp;[@ChangeDetection/](../MATLAB/+smi_stat/@ChangeDetection/README.md)              | change detection analysis methods
&nbsp;&nbsp;&nbsp;[@DiffusionEstimator/](../MATLAB/+smi_stat/@DiffusionEstimator/README.md)        | diffusion estimation methods
&nbsp;&nbsp;&nbsp;[@HMM/](../MATLAB/+smi_stat/@HMM/README.md)                   | hidden Markov model methods
[MATLAB/+smi_vis](../MATLAB/+smi_vis/README.md):
&nbsp;&nbsp;&nbsp;[@GenerateImages/](../MATLAB/+smi_vis/@GenerateImages/README.md)                 | general visualization of super-resolution data
&nbsp;&nbsp;&nbsp;[@GenerateMovies/](../MATLAB/+smi_vis/@GenerateMovies/README.md)                 | generate movies
&nbsp;&nbsp;&nbsp;[@InspectResults/](../MATLAB/+smi_vis/@InspectResults/README.md)                 | inspect super-resolution data
