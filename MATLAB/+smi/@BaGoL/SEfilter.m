function [Coords,Ind] = SEfilter(ROIs)
%SEfilter filters out NaNs and localizations with zero precisions 
% obj.NNDfilter(obj,SMD)
%
% INPUTS:
%   ROIs: structure containing X and Y coorinates. (nm)
%
% OUTPUT:
%   Coords: Structure containing filtered X and Y coordinates
%   Ind: Index of filtered localizations
%
% Created by:
%   Mohamadreza Fazel (Lidke Lab, 2022)
%
    [y,~]=size(ROIs.X);
     if y == 1
         XYT = [ROIs.X',ROIs.Y'];
         XYT_SE = [ROIs.X_SE',ROIs.Y_SE'];
     else
         XYT = [ROIs.X,ROIs.Y];
         XYT_SE = [ROIs.X_SE,ROIs.Y_SE];
     end  
        
     if isempty(XYT)
        Coords = [];
        Ind = [];
     else
         Ind1 = ~isnan(XYT(:,1)) & ~isnan(XYT(:,2));
         Ind2 = XYT_SE(:,1) ~=0 & XYT_SE(:,2) ~=0;
         Ind = Ind1 & Ind2;
         Coords = single(XYT(Ind,1:2));
     end
end