function [Coords,Ind] = NNDfilter(obj,ROIs)
%NNDfilter filters out localizations with less than a given number of neighbors
%   or more than a given number of neighbors.
% obj.NNDfilter(obj,SMD)
%
% NNDfilter takes a set of given coordinates and filters out localizations
% with less than a given number (NNN) or more than a given number (NNNmax)
% of neighbors within a given distance.
%
% INPUTS:
%   obj:    structure containing
%      NNR:    only accept at least NNN and no more than NNNmax localizations
%              within NNR (nm) for each localization
%      NNN:    Number of Nearest Neighbors
%      NNNmax: maximum Number of Nearest Neighbors
%   ROIs:   structure containing X and Y coordinates. (nm)
%
% OUTPUT:
%   Coords: Structure containing filtered X and Y coordinates
%   Ind:    Index of filtered localizations
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

     [m, ~] = size(XYT);

     if isempty(XYT)
        Coords = [];
        Ind = [];
     else
         if ~isempty(obj.NNN)
             % Are there sufficient points (m) to impose a maximum of NNNmax on
             % the number of nearest neighbors?
             if ~isempty(obj.NNNmax) && obj.NNNmax + 1 < m
                % Check self + NNs + 1 extra nearest neighbor beyond NNNmax.
                NNNn = max(obj.NNN, obj.NNNmax)+2;
                [~,Dis] = knnsearch(XYT,XYT,'k',NNNn);
                % NNN nearest neighbors (+ self) located within NNR.
                Ind = Dis(:,obj.NNN+1)<obj.NNR; 
                % NNNmax + self + 1 extra NN should NOT be located within NNR.
                Ind = Ind & Dis(:,obj.NNNmax+2)>obj.NNR;
                Coords = XYT(Ind,1:2);
             else
                % No checks needed on maximum number of nearest neighbors as
                % this restriction will be satisfied automatically because the
                % number of possibe neighbors is <= NNNmax.
                NNNn = obj.NNN+1;
                [~,Dis] = knnsearch(XYT,XYT,'k',NNNn);
                Ind = Dis(:,end)<obj.NNR; 
                Coords = XYT(Ind,1:2);
             end
         else
             Ind = ~isnan(XYT(:,1));
             Coords = single(XYT(:,1:2));
         end
     end
end
