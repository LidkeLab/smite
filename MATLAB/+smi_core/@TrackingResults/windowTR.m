function TRWindowed = windowTR(TR, MinFrameNum, MaxFrameNum, Verbose)
%windowTR abstracts a window of frame numbers for the output TR structure.
% Frames from MinFramNum to MaxFrameNum are retained and output to TRWindowed.
%
% INPUTS:
%    TR            TR structure
%    MinFrameNum   minimum frame number to select from TR
%    MaxFrameNum   maximum frame number to select from TR
%    Verbose       [default: false] if true, print nessage of what was done
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
   SubIndices = MinFrameNum <= SMD.FrameNum & SMD.FrameNum <= MaxFrameNum;
   SMDiso = smi_core.SingleMoleculeData.isolateSubSMD(SMD, SubIndices);
   if Verbose
      fprintf('Selected %d frames from %d frames\n',...
               numel(SMDiso.FrameNum), numel(SMD.FrameNum));
   end
   TRWindowed = smi_core.TrackingResults.convertSMDToTR(SMDiso);

end
