function assignROIs(obj,ROIs)
%assignROIs applies SE_Adjust and filters zero precisions or nans in the data   
%
% INPUTS:
%    ROIs: 2D Array of SMD structures (See genROIs())
%
% OUTPUT:
%    obj: BaGoL object populated with cluster coordinates in field ClusterSMD.
%      ClusterSMD: 
%       X:        Vector of X localization positions (nm)(Mx1)
%       Y:        Vector of Y localization positions (nm)(Mx1)
%       Z:        Vector of Z localization positions (nm)(Mx1)(Optional)
%       X_SE:     Vector of X localization standard errors (nm)(Mx1)
%       Y_SE:     Vector of Y localization standard errors (nm)(Mx1)
%       Z_SE:     Vector of Z localization standard errors (nm)(Mx1)(Optional)
%       FrameNum: Vector of localization frame numbers (Mx1)
%       ix:       Index of the ROI that a cluster comes from along the X axis
%       iy:       Index of the ROI that a cluster comes from along the Y axis
%

% Created by:
%    Mohamadreza Fazel (LidkeLab 2022)

obj.ClusterSMD = [];
[Ys,Xs] = size(ROIs);
for iy = 1:Ys
    for ix = 1:Xs
        
       [Coords,Ind] = BaGoL.SEfilter(ROIs(iy,ix));
       
        if size(Coords,1) <= 1
            continue; 
        end
        
        Len = length(obj.ClusterSMD);
        TX_SE = ROIs(iy,ix).X_SE(Ind);
        TY_SE = ROIs(iy,ix).Y_SE(Ind);
        if ~isempty(ROIs(iy,ix).FrameNum)
            TFrame = ROIs(iy,ix).FrameNum(Ind);
        end
        
        obj.ClusterSMD(Len+1).X = Coords(:,1);
        obj.ClusterSMD(Len+1).Y = Coords(:,2);
        obj.ClusterSMD(Len+1).X_SE = TX_SE(:)+obj.SE_Adjust;
        obj.ClusterSMD(Len+1).Y_SE = TY_SE(:)+obj.SE_Adjust;
        if isempty(ROIs(iy,ix).Z)
            obj.ClusterSMD(Len+1).Z = [];
            obj.ClusterSMD(Len+1).Z_SE = [];
        else
            obj.ClusterSMD(Len+1).Z = ROIs(iy,ix).Z(:);
            obj.ClusterSMD(Len+1).Z_SE = ROIs(iy,ix).Z_SE(:);
        end
        if ~isempty(ROIs(iy,ix).FrameNum)
            obj.ClusterSMD(Len+1).FrameNum = TFrame(:);
        else
            obj.ClusterSMD(Len+1).FrameNum = [];
        end
            obj.ClusterSMD(Len+1).ix = ix;
            obj.ClusterSMD(Len+1).iy = iy;
    end
end

end
