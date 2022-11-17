function precluster(obj,ROIs)
%precluster Finds the independent clusters for analysis with RJMCMC
% obj.precluster(ROIs)
% 
% precluster loops through the subregions, applies the NND filtering 
% (if applicable) and then uses hierarchical clustering to find independent
% clusters that will be analyzed individually by BaGoL_RJMCMC. The
% SMD structures for each structure is then saved inside the ClusterSMD
% structure, which is a property of BaGoL, see below.
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
%    Mohamadreza Fazel (LidkeLab 2020)

obj.ClusterSMD = [];
[Ys,Xs] = size(ROIs);
for iy = 1:Ys
    for ix = 1:Xs
        
    %[Coords,Ind] = obj.NNDfilter(ROIs(iy,ix));
    [Coords,Ind] = smi.BaGoL.SEfilter(ROIs(iy,ix));

        %clustering - give back array of SMD - one per cluster
        % There are three functions for hierarchical clustering that we have to 
        % use in turn. The first one is pdist which calculates the 
        % distance between every possible pair of the data points. This is the 
        % input to the next function which is linkage. The output of linkage
        % is a data structure defining a hierarchical clustering tree based on
        % the pdist results.
        if size(Coords,1) <= 1
            continue; 
        end
        try
            Dist = pdist(Coords,'euclid');
            TC = linkage(Dist,'single');
            IndClust = cluster(TC,'Cutoff',obj.Cutoff,'Criterion','distance');
        catch
            IndClust = ones(size(Dist(:,1)),'single');
        end
                    
        % cluster is the third function that we need to use. The input to this
        % function is the output of the previous one. This function has several
        % methods that we need to decide which one is suitable for our purpose. 
        % For example, the method 'distance' uses the data based on the
        % distance between the data points.
        NClust = max(IndClust); 
        MeanX = zeros(NClust,1);
        MeanY = zeros(NClust,1);
        for nn = 1:NClust
            MeanX(nn) = mean(Coords(IndClust==nn,1));
            MeanY(nn) = mean(Coords(IndClust==nn,2));
        end
        Len = length(obj.ClusterSMD);
        TX_SE = ROIs(iy,ix).X_SE(Ind);
        TY_SE = ROIs(iy,ix).Y_SE(Ind);
        if ~isempty(ROIs(iy,ix).FrameNum)
            TFrame = ROIs(iy,ix).FrameNum(Ind);
        end
        for nn = 1:max(IndClust)
            Indtrue = IndClust==nn;
            obj.ClusterSMD(Len+nn).X = Coords(Indtrue,1);
            obj.ClusterSMD(Len+nn).Y = Coords(Indtrue,2);
            obj.ClusterSMD(Len+nn).X_SE = TX_SE(Indtrue)+obj.SE_Adjust;
            obj.ClusterSMD(Len+nn).Y_SE = TY_SE(Indtrue)+obj.SE_Adjust;
            if isempty(ROIs(iy,ix).Z)
                obj.ClusterSMD(Len+nn).Z = [];
                obj.ClusterSMD(Len+nn).Z_SE = [];
            else
                obj.ClusterSMD(Len+nn).Z = ROIs(iy,ix).Z(Indtrue);
                obj.ClusterSMD(Len+nn).Z_SE = ROIs(iy,ix).Z_SE(Indtrue);
            end
                if ~isempty(ROIs(iy,ix).FrameNum)
                    obj.ClusterSMD(Len+nn).FrameNum = TFrame(Indtrue);
                else
                    obj.ClusterSMD(Len+nn).FrameNum = [];
                end
                obj.ClusterSMD(Len+nn).ix = ix;
                obj.ClusterSMD(Len+nn).iy = iy;
        end
    end
end

end
