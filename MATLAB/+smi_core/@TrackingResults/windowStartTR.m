function TRWindowed = windowStartTR(TR, MaxStartingFrameNum, Verbose)
%windowTR abstracts a window of frame numbers for the output TR structure.
% Trajectories starting at a frame number <= MaxStartingFrameNum are retained
% and output to TRWindowed.
%
% INPUTS:
%    TR           TR structure
%    MaxStartingFrameNum   maximum starting frame number sllowed for a
%                          trajectory
%    Verbose      [default: false] if true, print nessage of what was done
%
% OUTPUT:
%    TRWindowed   TR structure containing only trajectories starting at a
%                 frame number <= MaxStartingFrameNum

% Created by
%    Michael Wester (Lidke Lab, 2023)

   if ~exist('Verbose', 'var')
      Verbose = false;
   end

   % Convert TR structures back and forth between SMD structures so can
   % more easily restrict the frame numbers considered.
   SMD = smi_core.TrackingResults.convertTRToSMD(TR);
   SubIndices = zeros(size(SMD.FrameNum), 'logical');
   j = 1;
   for i = 1 : numel(TR)
      n = numel(TR(i).FrameNum);
      if TR(i).FrameNum(1) <= MaxStartingFrameNum
         SubIndices(j : j + n - 1) = true;
      end
      j = j + n;
   end
   SMDiso = smi_core.SingleMoleculeData.isolateSubSMD(SMD, SubIndices);
   if Verbose
      fprintf('Selected %d frames from %d frames\n',...
               numel(SMDiso.FrameNum), numel(SMD.FrameNum));
   end
   TRWindowed = smi_core.TrackingResults.convertSMDToTR(SMDiso);

end
