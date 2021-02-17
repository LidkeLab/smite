function plotNND_PDF(obj,SaveDir)
%plotNND_PDF makes and saves the NND PDF histogram of MAPN coordinates.
% obj.plotNND(SaveDir)
%
% Calculates the nearest neighbor distribution from the MAPN coordinates and
% creates a PDF histogram, also displaying the random nearest neighbor PDF that
% would occur for the same density of points.  This plot is scaled on the
% x-axis in two different ways: by 99% of the data and 99% of the random PDF.
% If the optional argument SaveDir is provided, the plots produced are saved to
% that directory as PNGs.
%
% INPUTS:
%    SaveDir: Save directory (Optional)
%
% OUTPUT:
%     None
%

% Created by:
%    Michael Wester and Mohamadreza Fazel (LidkeLab 2021 and 2019)

NBins=30;
[~,Dis]=knnsearch([obj.MAPN.X,obj.MAPN.Y],[obj.MAPN.X,obj.MAPN.Y],'k',2);
Dis = Dis(:,2);
% Scale the histogram by 99% of the data.
P = prctile(Dis,99);

% PDF for a random NN distribution.  Details at
%    https://en.wikipedia.org/wiki/Mean_inter-particle_distance
% density (#/nm^2)
%Area = (max(obj.MAPN.X)-min(obj.MAPN.X)) * (max(obj.MAPN.Y)-min(obj.MAPN.Y));
Area = obj.PImageSize^2;   % PImageSize = Posterior Image Size (nm)
maxD = max(Dis(Dis < P));
rho = numel(obj.MAPN.X) / Area;   % (#/nm^2)
a = 1 / sqrt(pi * rho);   % (nm)
r = 0 : maxD/1000 : maxD;   % (nm)
PDF = 2/a * (r/a) .* exp(-(r/a).^2);

figure
hold on
h = histogram(Dis(Dis<P),NBins);
h.Normalization = 'PDF';
plot(r, PDF, 'r-', 'LineWidth', 2);
xlabel('NND (nm)','FontSize',15);
ylabel('PDF','FontSize',15)
legend('data', 'random')
hold off

if nargin > 1
   print(gcf,fullfile(SaveDir,'NNDScaledData'),'-dpng'); 
   Dis = double(Dis);
   if size(Dis,2) > 1
      Dis = Dis';
   end
   save(fullfile(SaveDir,'NND.txt'),'Dis','-ascii')
end

% Scale the histogram by 99% of the random curve.
% Note: integrate(PDF(r), r, 0, inf) = 1 and
%       integrate(PDF(r), r, 0, R) = - exp(-(r/a)^2)|r=R + exp(-(r/a)^2)|r=0
%                                  = 1 - exp(-(R/a)^2) =>
%       f = 1 - exp(-(R/a)^2) => R = a*sqrt(-log(1 - f))   for f in [0, 1]
f = 0.99;
R = a * sqrt(-log(1 - f));
r = 0 : R/1000 : R;   % (nm)
PDF = 2/a * (r/a) .* exp(-(r/a).^2);

figure
hold on
h = histogram(Dis(Dis<R),NBins);
h.Normalization = 'PDF';
plot(r, PDF, 'r-', 'LineWidth', 2);
xlabel('NND (nm)','FontSize',15);
ylabel('PDF','FontSize',15)
legend('data', 'random')
hold off

if nargin > 1
   print(gcf,fullfile(SaveDir,'NNDScaledRandom'),'-dpng'); 
end

end
