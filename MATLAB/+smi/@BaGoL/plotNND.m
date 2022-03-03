function plotNND(obj,SaveDir)
%plotNND Makes and saves the NND histogram of MAPN coordinates
% obj.plotNND(SaveDir)
%
%  Calculates the nearest neighbor distribution from the MAPN coordinates 
%  and creates a histogram. If SaveDir input is given, it also saves the
%  plot to that directory as a PNG. 
%
% INPUTS:
%    SaveDir: Save directory (Optional)
%
% OUTPUT:
%     None
%

% Created by:
%    Mohamadreza Fazel (LidkeLab 2020)

NBins=30;
[~,Dis]=knnsearch([obj.MAPN.X,obj.MAPN.Y],[obj.MAPN.X,obj.MAPN.Y],'k',2);
Dis = Dis(:,2);
P = prctile(Dis,99);
figure
hist(Dis(Dis<P),NBins)
xlabel('NND(nm)','FontSize',15);ylabel('Frequency','FontSize',15)
if nargin > 1
   print(gcf,fullfile(SaveDir,'NND'),'-dpng'); 
   Dis = double(Dis);
   if size(Dis,2) > 1
      Dis = Dis';
   end
   save(fullfile(SaveDir,'NND.txt'),'Dis','-ascii')
end
end
