function FitFrame = fitsPerFrame(SMD, DatasetIndex)
%fitsPerFrame finds the fits per frame for the localizations in the given SMD
% structure.
%
% INPUTS:
%    SMD            SMD structure
%    DatasetIndex   dataset index for fits per frame for a single dataset 
%
% OUTPUTS:
%    FitFrame       fits per frame discovered for this/these dataset/s

   Nloc_frame = [];
   % Number of localizations per frame
   if ~exist('DatasetIndex', 'var')
       range = 1:max(SMD.DatasetNum);
   else
       range = DatasetIndex;
   end
   for jj=range
       for ii=1:SMD.NFrames
           idx=find(SMD.FrameNum==ii & SMD.DatasetNum==jj);
           Nloc_frame{jj}(ii)=length(idx);
       end
   end
   FitFrame = Nloc_frame{1};
   for ii=2:length(Nloc_frame)
       FitFrame=cat(2, FitFrame, Nloc_frame{ii});
   end

end
