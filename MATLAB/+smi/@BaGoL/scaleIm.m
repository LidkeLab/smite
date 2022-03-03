function ImOut=scaleIm(Im,Percentile)
%scaleIm Scales and clips the image intensity to improve image contrast
% ImOut=BaGoL.scaleIm(Im,Percentile)
%
% The input image is linearly scaled between zero and a percentile of 
% intensity of non-zero pixels. Values above that percentile are set to the 
% value of the intensity at that percentile.  
%
% INPUT:
%    Im:         The input image
%    Percentile: Percentile used for scaling. (Default = 99)
%
% OUTPUT:
%    Imout:      Clipped and Scaled Image
%

% Created by:
%    Mohmadreza Fazel (LidkeLab 2019)
%
    if nargin < 2
        Percentile = 99;
    end
    ImOut=Im;
    P=prctile(Im(Im>0),Percentile);
    ImOut(ImOut>P)=P;
    ImOut=255*ImOut/P;
end
