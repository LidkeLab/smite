function [Report]=optimPSFZernike_unitTest()  
%optimPSFZernike_unitTest tests optimPSFZernike functionality.
%
% OUTPUTS:
%   Report
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%
 
%% Blind start. Search for Astigmatism
P=smi_psf.PointSpreadFunction.createPSFStruct();
PhaseMask=0;
PhaseMask(6)=1; %Vertical Astigmatism 
clc;
[PFound,PSF,CRLB_Astig]=smi_psf.PointSpreadFunction.optimPSFZernike(P,PhaseMask);

%%  Directed Start. Search for Tetrapod
P=smi_psf.PointSpreadFunction.createPSFStruct();
PhaseMask=0;
PhaseMask(6)=1; %Vertical Astigmatism 
PhaseMask(12)=1; %Vertical Second Order Astigmatism
PhaseMask(24)=0; %Vertical Second Order Astigmatism
P.ZC_Phase(6)=1;
P.ZC_Phase(12)=-2;
P.ZC_Phase(24)=2;
clc;
[PFound_Tet,PSF3,CRLB_Tetra]=smi_psf.PointSpreadFunction.optimPSFZernike(P,PhaseMask,P.ZC_Phase);
PFound_Tet.ZC_Phase

end
