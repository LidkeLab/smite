function [Report] = psfROIStack_unitTest()
%psfROIStack_unitTest Tests psfROIStack functionality.

Report = 0;

SMD.X = rand([1,100000])*5+5;
SMD.Y = rand([1,100000])*5+5;
SMD.Z = rand([1,100000])*10-5;
SMD.Photons = 1000*ones([1,100000]);
SMD.Bg = 15*ones([1,100000]);
SZ = 15;
XYSamPerPix = 6;
ZSamPerUnit = 20;
Z0=1;           
d=1;
PSF0=0.7;

[XGr,YGr,ZGr]=meshgrid(-20+1/(2*XYSamPerPix):1/XYSamPerPix:20-1/(XYSamPerPix*2),...
               -20+1/(XYSamPerPix*2):1/XYSamPerPix:20-1/(2*XYSamPerPix),-15:1/ZSamPerUnit:15);

PSF = (normpdf(XGr,0,PSF0*sqrt(1+((Z0-ZGr)/d).^2)).*...
       normpdf(YGr,0,PSF0*sqrt(1+((Z0+ZGr)/d).^2)));
tic;
[Model,Data]=smi_psf.PointSpreadFunction.psfROIStack(PSF,XYSamPerPix,ZSamPerUnit,SZ,SMD);
T = toc;
fprintf('psfROIStack is successfully tested.\n');
fprintf('%g blobs of 3D super-resolution data were generated in %g seconds.\n',size(SMD.X,2),T);

Report = 1;
end
