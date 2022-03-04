function Ind = removeOverlap(obj,ROIs,X,Y,ii)
%removeOverlap Finds coordinates located in overlapping regions.
% Ind = obj.removeOverlap(ROIs,X,Y,ii)
%
% removeOverlap() finds the coordinates in the given list that are inside 
% the overlapping regions between the adjacent ROIs and returns the index 
% of the coordinates that are not inside the overlapping regions. obj is 
% not modified by this method.
%
% Uses the class propoerty 'Overlap' and the location of the ROI to which
% the cluster of index ii belongs to make the calculation. 
%
% INPUTS:
%    ROIs:  2D Array of SMD structures (See genROIs())
%    X:     Vector of X coordinates (Mx1)
%    Y:     Vector of Y coordinates (Mx1)
%    ii:    Index of cluster 
%
% OUTPUT:
%    Ind:   The indices of the emitters that are inside the non-overlapping
%           regions.
%

% Created by:
%    Mohamadreza Fazel (LidkeLab 2019)

ix = obj.ClusterSMD(ii).ix;
iy = obj.ClusterSMD(ii).iy;
XPix = max(obj.SMD.X);
YPix = max(obj.SMD.Y);
[Ysub,Xsub] = size(ROIs); 
MinX = min(obj.SMD.X);
MinY = min(obj.SMD.Y);
if MinX < obj.ROIsize
    MinX = 0;
end
if MinY < obj.ROIsize
    MinY = 0; 
end
sy = obj.ROIsize*(iy-1);
ey = obj.ROIsize*iy;
sx = obj.ROIsize*(ix-1);
ex = obj.ROIsize*ix;
if iy == Ysub
    ey = YPix;
end
if ix == Xsub
    ex = XPix;
end                    
Ind = find(X - MinX >= sx & X - MinX < ex ...
      & Y - MinY >= sy & Y - MinY < ey);
end
