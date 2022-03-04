function errPlot(SMD) 
%errPlot Plots the input precisions, where circles represent the errors.
% BaGoL.errPlot(SMD) 
%
% The plot contains circles that represent localizations within the SMD. 
% For each localization, the center of the circle is located at X,Y and 
% the radii of circles are 2 sigma where sigma is taken as
% sigma=sqrt(X_SE^2+Y_SE^2). 
%
% INPUT:
%    SMD: Structure containing the following fields:
%       X:      Vector of X localization positions (nm)(Nx1)
%       Y:      Vector of Y localization positions (nm)(Nx1)
%       X_SE:   Vector of X localization standard errors (nm)(Nx1)
%       Y_SE:   Vector of Y localization standard errors (nm)(Nx1)
%
% OUTPUT:
%    None.
%

% Created by:
%    Mohamadreza Fazel (LidkeLab 2019)

Ang = 0:0.1:2*pi+0.25;
XCir = zeros(length(Ang),length(SMD.X),'single');
YCir = XCir;
for nn = 1:length(SMD.X)
    Prec = 2*sqrt((SMD.X_SE(nn)^2+SMD.Y_SE(nn)^2)/2);
    XCir(:,nn) = Prec*cos(Ang')+SMD.X(nn);
    YCir(:,nn) = Prec*sin(Ang')+SMD.Y(nn);
end
figure;
plot(XCir,YCir,'m','linewidth',1.1)
axis equal
end
