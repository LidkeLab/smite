function plotMAPN(obj,PlotType,SaveDir,PlotVis)
%plotCluster: Plots the MAPN positions together with input localizations
% obj.plotCluster(FlagOrDir,FlagPlot)
%
% The first plot is color-coded frame numbers and the MAPN position.
% The second plot is plot of data points before and after filtering and the
% MAPN coordinates.  The plots are saved in the directory given (optional)
% in fig-format.
%
% INPUTS:
%   PlotType: Type of plot to generate (Default=3)  
%           1: Clusters color coded with MAPN locations
%           2: Shows filter and unfiltered localizations with MAPN
%           3: Both 1 and 2
%   SaveDir: Save directory (Optional) 
%   PlotVis: Show the plots (Default=off) (Optional)
%
% OUTPUT:
%    None
%

% Created by:
%    Mohamadreza Fazel (LidkeLab 2019)
%

if nargin < 4
   PlotVis = 'off'; 
end
if nargin < 2
   PlotType = 3; 
end

if PlotType == 1 || PlotType == 3
    figure('Visible',PlotVis);hold;
    P1 = scatter(obj.ClusterSMD(1).X,obj.ClusterSMD(1).Y,[],obj.ClusterSMD(1).FrameNum);
    for mm = 2:length(obj.ClusterSMD)
        scatter(obj.ClusterSMD(mm).X,obj.ClusterSMD(mm).Y,[],obj.ClusterSMD(mm).FrameNum)
    end
    P2 = plot(obj.MAPN.X,obj.MAPN.Y,'*k');
    colorbar();xlabel('X(nm)','FontSize',16);ylabel('Y(nm)','FontSize',16)
    axis equal;
    legend([P1,P2],{'Localizations','Groups centers'})    
    if nargin > 2
       saveas(gcf,fullfile(SaveDir,'LocsScatter-MAPN.fig')) 
    end
end
if PlotType == 2 || PlotType ==3
    figure('Visible',PlotVis);hold;
    P1 = plot(obj.SMD.X,obj.SMD.Y,'.k');
    P2 = plot(obj.ClusterSMD(1).X,obj.ClusterSMD(1).Y,'.','color',[0,0.447,0.741]);
    for mm = 2:length(obj.ClusterSMD)
        plot(obj.ClusterSMD(mm).X,obj.ClusterSMD(mm).Y,'.','color',[0,0.447,0.741]) 
    end
    P3 = plot(obj.MAPN.X,obj.MAPN.Y,'*','color',[0.3,0.15,0.15],'linewidth',1.25);
    axis equal;xlabel('X(nm)','FontSize',16);ylabel('Y(nm)','FontSize',16)
    legend([P1,P2,P3],{'Before filtering','After filtering','Groups centers'})
    if nargin > 2
       saveas(gcf,fullfile(SaveDir,'Signal-BG-Cluster.fig')) 
    end
end
end
