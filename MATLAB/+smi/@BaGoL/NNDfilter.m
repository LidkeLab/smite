function [Coords,Ind] = NNDfilter(obj,ROIs)
%NNDfilter filters out localizations with less than a given number of neighbors
% obj.NNDfilter(obj,SMD)
%
% NNDfiler takes a set of given coordinates and filters out localizations
% with less than a given number of neighbors within a given distance.
%
% INPUTS:
%   ROIs: structure containing X and Y coorinates. (nm)
%
% OUTPUT:
%   Coords: Structure containing filtered X and Y coordinates
%   Ind: Index of filtered localizations
%
% Created by"
%   Mohamadreza Fazel (Lidke Lab, 2020)
%
    [y,~]=size(ROIs.X);
     if y == 1
         XYT = [ROIs.X',ROIs.Y'];
     else
         XYT = [ROIs.X,ROIs.Y];
     end  
        
     if isempty(XYT)
        Coords = [];
        Ind = [];
     else
         if ~isempty(obj.NNN)
             NNNn = obj.NNN+1;
             [~,Dis] = knnsearch(XYT,XYT,'k',NNNn);
             Ind = Dis(:,end)<obj.NNR; 
             Coords = XYT(Ind,1:2);
         else
             Ind = ~isnan(XYT(:,1));
             Coords = single(XYT(:,1:2));
         end
     end
end