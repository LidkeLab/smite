function TRWindowed = windowTR(TR, MinFrameNum, MaxFrameNum, Verbose)
%windowTR abstracts a window of frame numbers for the output TR structure.
% Frames from MinFramNum to MaxFrameNum are retained and output to TRWindowed.
%
% INPUTS:
%    TR            TR structure
%    MinFrameNum   minimum frame number to select from TR
%                  if < 0, it is the -minFrameNum of frames to select, e.g.,
%                  -3 means select the first 3 frames of each trajectory
%    MaxFrameNum   maximum frame number to select from TR
%                  if > 0 and < 1, it is the fraction of total frames to
%                  select, e.g. 0.3 means select the first 30% of the frames of
%                  each trajectory
%    Verbose       [default: false] if true, print nessage of what was done
% NOTE: if MinFrameNum < 0 AND MaxFrameNum is in (0, 1), MinFrameNum action
%       takes precedence.
%
% OUTPUT:
%    TRWindowed    TR structure containing only frames numbered between
%                  MinFrameNum and MaxFrameNum

% Created by
%    Michael Wester (Lidke Lab, 2023)

   if ~exist('Verbose', 'var')
      Verbose = false;
   end

   % Convert TR structures back and forth between SMD structures so can
   % more easily restrict the frame numbers considered.
   SMD = smi_core.TrackingResults.convertTRToSMD(TR);
   if MinFrameNum >= 1 && MaxFrameNum >= 1
      SubIndices = MinFrameNum <= SMD.FrameNum & SMD.FrameNum <= MaxFrameNum;
   else
      % NTraj is the total number of trajectories in the TR structure.
      NTraj = numel(TR);
      % LTraj counts the number of frames per trajectory in the TR structure.
      LTraj = arrayfun(@(i) numel(TR(i).FrameNum), 1:NTraj, ...
                       'UniformOutput', false);
      LTraj = arrayfun(@(i) LTraj{i}, 1:NTraj);
      NFramesTotal = sum(LTraj);
      SubIndices = zeros(1, NFramesTotal, 'logical');

      if MinFrameNum < 0
         % Select an absolute number of frames for each trajectory (ifi
         % available).
         MinFrameCount = -MinFrameNum;
         MTraj = arrayfun(@(i) min(LTraj(i), MinFrameCount), 1:NTraj);
      elseif MaxFrameNum > 0 && MaxFrameNum < 1
         % Select a fraction of the total frames in the TR structure for each
         % trajectory starting at the first frame.  Make sure that at least one
         % frame is selected.
         MaxFraction = MaxFrameNum;
         MTraj = max(1, round(MaxFraction .* LTraj));
      end

      for i = 1 : NTraj
         j = sum(LTraj(1 : i - 1)) + 1;
         SubIndices(j : j + MTraj(i) - 1) = true;
      end   
   end
   SMDiso = smi_core.SingleMoleculeData.isolateSubSMD(SMD, SubIndices);
   if Verbose
      fprintf('Selected %d frames from %d frames\n',...
               numel(SMDiso.FrameNum), numel(SMD.FrameNum));
   end
   TRWindowed = smi_core.TrackingResults.convertSMDToTR(SMDiso);

end
