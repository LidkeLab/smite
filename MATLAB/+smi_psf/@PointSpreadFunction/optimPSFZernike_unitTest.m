function [Report]=optimPSFZernike_unitTest()  
%optimPSFZernike_unitTest Tests optimPSFZernike functionality.
%
% OUTPUTS:
%   Report
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%

Report = 0;

SaveDir = smi_helpers.mkSMITETmpDir('unitTest', 'optimPSFZernike');
 
%% Blind start. Search for Astigmatism
P=smi_psf.PointSpreadFunction.createPSFStruct();
PhaseMask=0;
PhaseMask(6)=1; %Vertical Astigmatism 
%clc;
[PFound,PSF,CRLB_Astig]=smi_psf.PointSpreadFunction.optimPSFZernike(P,PhaseMask);
saveas(gcf, fullfile(SaveDir, 'oPZ1.png'));

%%  Directed Start. Search for Tetrapod
P=smi_psf.PointSpreadFunction.createPSFStruct();
PhaseMask=0;
PhaseMask(6)=1; %Vertical Astigmatism 
PhaseMask(12)=1; %Vertical Second Order Astigmatism
PhaseMask(24)=0; %Vertical Second Order Astigmatism
P.ZC_Phase(6)=1;
P.ZC_Phase(12)=-2;
P.ZC_Phase(24)=2;
%clc;
[PFound_Tet,PSF3,CRLB_Tetra]=smi_psf.PointSpreadFunction.optimPSFZernike(P,PhaseMask,P.ZC_Phase);
PFound_Tet.ZC_Phase
saveas(gcf, fullfile(SaveDir, 'oPZ2.png'));

Report = 1;

end
