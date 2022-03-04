function plotMAPN(obj,SaveDir,PlotVis)
%plotCluster: Plots the MAPN positions together with input localizations
% obj.plotCluster(FlagOrDir,FlagPlot)
%
% The first plot is color-coded frame numbers and the MAPN position.
% The second plot is plot of data points before and after filtering and the
% MAPN coordinates.  The plots are saved in the directory given (optional)
% in fig-format.
%
% INPUTS:
%   SaveDir: Save directory (Optional) 
%   PlotVis: Show the plots (Default=off) (Optional)
%
% OUTPUT:
%    None
%

% Created by:
%    Mohamadreza Fazel (LidkeLab 2019)
%

if nargin < 3
   PlotVis = 'off'; 
end

    figure('Visible',PlotVis);hold;
    P1 = scatter(obj.ClusterSMD(1).X,obj.ClusterSMD(1).Y,[],obj.ClusterSMD(1).FrameNum);
    for mm = 2:length(obj.ClusterSMD)
        scatter(obj.ClusterSMD(mm).X,obj.ClusterSMD(mm).Y,[],obj.ClusterSMD(mm).FrameNum)
    end
    P2 = plot(obj.MAPN.X,obj.MAPN.Y,'*k');
    colorbar();xlabel('X(nm)','FontSize',16);ylabel('Y(nm)','FontSize',16)
    axis equal;
    legend([P1,P2],{'Localizations','Groups centers'})    
    saveas(gcf,fullfile(SaveDir,'LocsScatter-MAPN.fig')) 

end
