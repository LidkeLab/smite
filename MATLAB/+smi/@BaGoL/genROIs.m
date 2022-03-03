function ROIs = genROIs(obj)
%genROIs takes the input coordinates and splits them into smaller regions.
% ROIs = obj.genROIs()
%
% genROIs() uses the ROIsize and overlap size (see BaGoL properties) to
% break down the input localizations into smaller regions of interest.
% These regions are overlapped with their adjacent regions to avoid edge
% artifacts. The coordinates, precisions, frame numbers associated to each
% ROI are saved inside an element of the output structure array.
%
% INPUT:
%       X:        Vector of X localization positions (nm)(Nx1)
%       Y:        Vector of Y localization positions (nm)(Nx1)
%       Z:        Vector of Z localization positions (nm)(Nx1)(Optional)
%       X_SE:     Vector of X localization standard errors (nm)(Nx1)
%       Y_SE:     Vector of Y localization standard errors (nm)(Nx1)
%       Z_SE:     Vector of Z localization standard errors (nm)(Nx1)(Optional)
%       FrameNum: Vector of localization frame numbers (Nx1)
%
% OUTPUT:
%    ROIs: 2D Array of SMD structures with fields:
%       X:        Vector of X localization positions (nm)(Mx1)
%       Y:        Vector of Y localization positions (nm)(Mx1)
%       Z:        Vector of Z localization positions (nm)(Mx1)(Optional)
%       X_SE:     Vector of X localization standard errors (nm)(Mx1)
%       Y_SE:     Vector of Y localization standard errors (nm)(Mx1)
%       Z_SE:     Vector of Z localization standard errors (nm)(Mx1)(Optional)
%       FrameNum: Vector of localization frame numbers (Mx1)
%       Xsize:    Size of the region along X (nm)  
%       Ysize:    Size of the region along Y (nm)
%

% Created by:
%    Mohamadreza Fazel (LidkeLab 2019)

XPix = max(obj.SMD.X);
YPix = max(obj.SMD.Y);
MinX = min(obj.SMD.X);
MinY = min(obj.SMD.Y);
if MinX < obj.ROIsize && MinX > 0
    MinX = 0;
end
if MinY < obj.ROIsize && MinY > 0
    MinY = 0; 
end
Xsub = ceil((XPix-MinX)/obj.ROIsize); %number of the sub-regions along the X-axis
Ysub = ceil((YPix-MinY)/obj.ROIsize); %number of the sub-regions along the Y-axis
ROIs = [];
for ix = 1:Xsub
    for iy = 1:Ysub
        ROIex = obj.ROIsize + 2*obj.Overlap;
        ROIey = obj.ROIsize + 2*obj.Overlap;
        sy = obj.ROIsize*(iy-1)-obj.Overlap;
        ey = obj.ROIsize*iy+obj.Overlap;
        sx = obj.ROIsize*(ix-1)-obj.Overlap;
        ex = obj.ROIsize*ix+obj.Overlap;
        if iy == 1 
            ROIey = obj.ROIsize + obj.Overlap;
            sy = 0;
            ey = obj.ROIsize + obj.Overlap;
        end
        if ix == 1
            ROIex = obj.ROIsize + obj.Overlap;
            sx = 0;
            ex = obj.ROIsize + obj.Overlap;
        end
        if iy == Ysub
            %ey = YPix;
            ROIey = YPix-obj.ROIsize*(Ysub-1)+obj.Overlap;
        end
        if ix == Xsub
            %ex = XPix;
            ROIex = XPix-obj.ROIsize*(Xsub-1)+obj.Overlap;
        end
        if Xsub == 1
            ROIex = XPix;
        end
        if Ysub == 1
            ROIey = YPix;
        end
        Ind = obj.SMD.X-MinX >= sx & obj.SMD.X-MinX <= ex ...
             & obj.SMD.Y-MinY >= sy & obj.SMD.Y-MinY <= ey;
        ROIs(iy,ix).X = obj.SMD.X(Ind);
        ROIs(iy,ix).Y = obj.SMD.Y(Ind);
        
        if isfield(obj.SMD,'Z') 
            if ~isempty(obj.SMD.Z)
                ROIs(iy,ix).Z = obj.SMD.Z(Ind);
            else
                ROIs(iy,ix).Z = [];
            end
        else
            ROIs(iy,ix).Z = [];
            obj.SMD.Z=[];
        end
        
        ROIs(iy,ix).X_SE = obj.SMD.X_SE(Ind);
        ROIs(iy,ix).Y_SE = obj.SMD.Y_SE(Ind);
        
        if isfield(obj.SMD,'Z') 
            if ~isempty(obj.SMD.Z)
                ROIs(iy,ix).Z_SE = obj.SMD.Z_SE(Ind);
            else
                ROIs(iy,ix).Z_SE = [];
            end
        else
            ROIs(iy,ix).Z_SE = [];
        end
        
        if ~isempty(obj.SMD.FrameNum)
            ROIs(iy,ix).FrameNum = obj.SMD.FrameNum(Ind);
        else
            ROIs(iy,ix).FrameNum = []; 
        end
        ROIs(iy,ix).Xsize = ROIex;
        ROIs(iy,ix).Ysize = ROIey;
    end
end
end
