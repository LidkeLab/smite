function FitsFrame = fitsPerFrame(SMD)
%fitsPerFrame finds the fits per frame for the localizations in the given SMD
% structure.

   Nloc_frame = [];
   % Number of localizations per frame
   for jj=1:max(SMD.DatasetNum)
       for ii=1:max(SMD.FrameNum)
           idx=find(SMD.FrameNum==ii & SMD.DatasetNum==jj);
           Nloc_frame{jj}(ii)=length(idx);
       end
   end
   FitFrame = Nloc_frame{1};
   for ii=1:length(Nloc_frame)-1
       FitFrame=cat(2, FitFrame, Nloc_frame{ii+1});
   end

end
