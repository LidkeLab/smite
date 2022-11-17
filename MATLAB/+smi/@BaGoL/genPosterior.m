function [PostIm] = genPosterior(obj,PostIm,SZ,Chain,ROIs,ii)
%genPosterior Updates Posterior Image using RJMCMC chain for a subregion (ROI)
% [PostIm] = obj.genPosterior(PostIm,SZ,Chain,ROIs,ii)
%
% genPosterior() updates the Posterior Image with emitters from all the 
% states of the given chain. The coordinates that overlap with other ROIs 
% (see genROIs()) are first eliminated from the list of coordinates. The 
% final list of coordinates are then used to update the posterior image. 
% For each coordinate in the chain, the corresponding pixel in the
% Posterior Image is incremented by one. Therefore, the Posterior Image is
% the histogram image of the chain.
%
% INPUTS:
%    PostIm: Posterior Image to be updated
%    SZ:     Size of Posterior Image (nm)[Scalar]
%    Chain:  RJMCMC chain for the cluster (See BaGoL_RJMCMC())
%    ROIs:   2D Array of SMD structures (see genROIs())
%    ii:     Index of cluster associated with the input chain.
%
% OUTPUTS: 
%    PostIm: Updated Posterior Image
%

% Created by:
%    Mohamadreza Fazel (LidkeLab 2019)

XChain = [];
YChain = [];

for nn = 1:length(Chain)
    if ~isempty(Chain(nn).X)
        XChain = cat(1,XChain,Chain(nn).X(:));
        YChain = cat(1,YChain,Chain(nn).Y(:));
    end
end
Ind = obj.removeOverlap(ROIs,XChain,YChain,ii);

if isempty(obj.XStart)
    MinX = -((obj.PImageSize - (max(obj.SMD.X)-min(obj.SMD.X)))/2 - min(obj.SMD.X)); 
    MinY = -((obj.PImageSize - (max(obj.SMD.Y)-min(obj.SMD.Y)))/2 - min(obj.SMD.Y));
else
    MinX = -obj.XStart;
    MinY = -obj.YStart;
end
YInd = floor((YChain(Ind)+MinY)/obj.PixelSize);
XInd = floor((XChain(Ind)+MinX)/obj.PixelSize);

IndRemove = XInd>SZ | XInd < 1 | YInd>SZ | YInd < 1;
XInd(IndRemove) = [];
YInd(IndRemove) = [];

LinInd = sub2ind([SZ,SZ],YInd,XInd);
[Counts,Pixel]=hist(LinInd,unique(LinInd));
IndRemove = Counts == 0;
Counts(IndRemove) = [];
Pixel(IndRemove) = [];
PostIm(Pixel) = PostIm(Pixel) + Counts';

end
