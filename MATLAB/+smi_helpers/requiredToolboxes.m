function requiredToolboxes()
% Print out required toolboxes for each directory in the SMITE directory
% structure.

% Created by
%    Michael J. Wester (2021, LidkeLab)

SMITEdirs = {
   'ExternalSoftware'
   'ExternalSoftware/FRCresolution_software'
   'ExternalSoftware/PlotSpread'
   'ExternalSoftware/uipickfiles'
   'MATLAB'
   'MATLAB/+smi'
   'MATLAB/+smi/@BaGoL'
   'MATLAB/+smi/@Publish'
   'MATLAB/+smi/@SMLM'
   'MATLAB/+smi/@SPT'
   'MATLAB/+smi_cluster'
   'MATLAB/+smi_cluster/@Clustering'
   'MATLAB/+smi_cluster/@PairCorrelation'
   'MATLAB/+smi_cluster/@StatisticsClustering'
   'MATLAB/+smi_core'
   'MATLAB/+smi_core/@ChannelRegistration'
   'MATLAB/+smi_core/@DataToPhotons'
   'MATLAB/+smi_core/@DriftCorrection'
   'MATLAB/+smi_core/@FRC'
   'MATLAB/+smi_core/@FindROI'
   'MATLAB/+smi_core/@FrameConnection'
   'MATLAB/+smi_core/@LoadData'
   'MATLAB/+smi_core/@LocalizeData'
   'MATLAB/+smi_core/@SingleMoleculeData'
   'MATLAB/+smi_core/@SingleMoleculeFitting'
   'MATLAB/+smi_core/@Threshold'
   'MATLAB/+smi_core/@TrackingResults'
   'MATLAB/+smi_helpers'
   'MATLAB/+smi_helpers/@ROITools'
   'MATLAB/+smi_psf'
   'MATLAB/+smi_psf/@PointSpreadFunction'
   'MATLAB/+smi_psf/@Zernike'
   'MATLAB/+smi_sim'
   'MATLAB/+smi_sim/@GaussBlobs'
   'MATLAB/+smi_sim/@SimSMLM'
   'MATLAB/+smi_sim/@SimSPT'
   'MATLAB/+smi_stat'
   'MATLAB/+smi_stat/@DiffusionEstimator'
   'MATLAB/+smi_stat/@HMM'
   'MATLAB/+smi_vis'
   'MATLAB/+smi_vis/@GenerateMovies'
   'MATLAB/+smi_vis/@InspectResults'
   'MATLAB/+smi_vis/@GenerateImages'
   'MATLAB/examples'
   'MATLAB/mex'
   'MATLAB/ptx'
   'MATLAB/source'
   'MATLAB/source/c'
   'MATLAB/source/cuda'
   'MATLAB/source/cuda/smi_cuda_FindROI'
   'MATLAB/source/cuda/smi_cuda_PSFSample3DBlob'
   'MATLAB/source/cuda/smi_cuda_gaussBlobROIStack'
   'MATLAB/source/cuda/smi_cuda_gaussMLEv2'
};

nSMITEdirs = numel(SMITEdirs);
for i = 1 : nSMITEdirs
   SMITEdir = fullfile(userpath, 'smite', SMITEdirs{i});
   [fList, pList] = matlab.codetools.requiredFilesAndProducts(SMITEdir);
   nProducts = numel(pList);
   fprintf('%s:\n', SMITEdirs{i});
   for j = 1 : nProducts
      % Eliminate non-toolboxes from the output.
      if ~isempty(strfind(pList(j).Name, ' Toolbox'))
         fprintf('   %s\n', pList(j).Name);
      end
   end
end

end
