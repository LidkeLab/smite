function [Model, Data] = genImageStack(obj)
%genImageStack generates image stacks without and with Poisson noise.
%
% INPUTS:
%    obj            SimSMLM object
%       SMD_Model   labeled true coordinates with blinks
%       PSFSigma    point spread function sigma size (Pixels)
%       Bg          background count rate (counts/pixel)
%
% OUTPUTS:
%    Model          Gaussian blob image stack that is noiseless
%    Data           [OPTIONAL] Gaussian blob image stack corrupted with Poisson
%                   noise

% Created by:
%    Sajjad Khan and Michael Wester (2021, LidkeLab)

   % Generate the blobs without Poisson noise.  Copy the SMD structure.
   SMD_Model = obj.SMD_Model;
   % Below needed for generating blob images.
   SMF = smi_core.SingleMoleculeFitting();
   SMF.BoxFinding.BoxSize = ceil(4 * 2 * obj.PSFSigma);
   %SMD_Model.FitBoxSize = ceil(4 * 2 * obj.PSFSigma);
   % Temporarily convert FrameNum into an absolute frame number for the call
   % to gaussBlobImage.  Need a less resource intensive way to incorporate
   % NDatasets > 1.
   NFrames = obj.NDatasets * obj.NFrames;
   %FrameNum = SMD_Model.FrameNum;
   SMD_Model.NDatasets = 1;
   SMD_Model.NFrames = NFrames;
   SMD_Model.FrameNum = (SMD_Model.DatasetNum - 1) * obj.NFrames + FrameNum;
   [Model] = smi_sim.GaussBlobs.gaussBlobImage(SMD_Model, SMF);
   %[Model] = smi_sim.GaussBlobs.gaussBlobImage(obj.SZ, NFrames, SMD_Model, ...
   %                                            0, 0, 0);
   % Add in background.
   Model = Model + obj.Bg;

   if nargout > 1
      % Corrupt with Poisson noise.
      Data = poissrnd(single(Model));
      NoisyImage = zeros(size(Data(:, :, 1)));
      Data = Data + randn(size(Data)) .* repmat(NoisyImage, [1, 1, NFrames]);
   end

end
