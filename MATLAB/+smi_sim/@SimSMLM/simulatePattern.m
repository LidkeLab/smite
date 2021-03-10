function [SMD_True, SMD_True_Labeled, SMD_Model, SMD_Model_Noisy, Model, ...
          Data] = simulatePattern(obj, pattern)

   switch pattern
   case 'SiemensStar'
      % Siemen's star
      if isempty(obj.NWings)
         error('Siemen''s star must define NWings!');
      end
      SMD_True = obj.simStar(obj.NWings);
   case 'kTets'
      % Generate k-tets in the simulation region (units are pixels).
      if isempty(obj.OrderkTet) || isempty(obj.RadiuskTet)
         error('kTets must define OrderkTet and RadiuskTet!');
      end
      SMD_True = obj.kTets(obj.OrderkTet, obj.RadiuskTet);
   otherwise
      error('Unknown pattern!');
   end
   if nargout == 1 return; end

   % Apply labeling efficiency.
   SMD_True_Labeled = obj.applyLabelEffic(SMD_True);
   if nargout <= 2 return; end

   % Generate blinks (units are pixels).
   SMD_Model = obj.genBlinks(SMD_True_Labeled, obj.StartState);
   if nargout <= 3 return; end

   % EITHER, generate an SMD structure with positional and intensity noise.
   SMD_Model_Noisy = obj.genNoisySMD(SMD_Model);
   if nargout <= 4 return; end

   % OR ALTERNATIVELY, generate the blobs without Poisson noise (units are
   % pixels).  Below needed for generating blob images.
   SMD_Model.FitBoxSize = ceil(4 * 2 * obj.PSFSigma);
   % Temporarily convert FrameNum into an absolute frame number for the call
   % to gaussBlobImage.
   FrameNum = SMD_Model.FrameNum;
   SMD_Model.FrameNum = (SMD_Model.DatasetNum - 1) * obj.NFrames + FrameNum;
   Model = smi_sim.GaussBlobs.gaussBlobImage(obj.SZ,                      ...
                                             obj.NDatasets * obj.NFrames, ...
                                             SMD_Model, 0, 0, 0);
   SMD_Model.FrameNum = FrameNum;
   if nargout <= 5 return; end

   % Generate the blobs having Poisson noise.
   Data = obj.genNoisyData(Model);
end
