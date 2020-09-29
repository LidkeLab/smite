function [PSFStruct,PSF,CRLB]=optimPSFZernike(PSFStruct,PhaseMask,StartPhase,Photons,Bg)  
%optimPSFZernike Optimize PSF via Zernike Coef search. 
%   The PSF is generated using the following model:
%   PSF = |F[OTF]|^2
%   % more...
%
%   Noll ordering is used as a linear index.    
%
%   PSF is normalized such that the integral over all space = 1
%
% INPUTS:
%   PSFStruct:  PSF Structure.  (Default = createPSFStruct())
%   PhaseMask:  Nx1 vector of indicating what Zernike Coef will be adjusted
%   StartPhase: Nx1 vector of Coef as starting positions (Default = blind search)    
%   Photons:    Mx1 vector where M is length(P.Z)
%   Bg:         Mx1 vector where M is length(P.Z)
%
% OUTPUTS:
%   PSFStruct:  Updated PSF Structure.  
%   PSF:        Found PSF
%   CRLB:       CRLB at each z position
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%
    
if nargin <1
    P=smi_psf.PointSpreadFunction.createPSFStruct();
else
    P=PSFStruct;
end

if nargin <2
    PhaseMask=[0;0;0;ones(12,1)]
end

NCoef=sum(PhaseMask);
ZC0=ones(NCoef,1);

if nargin<4
   Photons=1000;
   Bg=10;
end

if nargin<3 %Blind search
    FValOld=costCRLB(ZC0,P,Photons,Bg,PhaseMask);
    ZCF_Best=ZC0;
    NAttempts=50;
    
    for nn=1:NAttempts
        ZC0=2*randn(NCoef,1); %first 3 are piston, and tilt, which don't help.
        Fval= costCRLB(ZC0,P,Photons,Bg,PhaseMask);
        if Fval<FValOld
            ZCF_Best=ZC0;
            FValOld=Fval;
        end
    end
else
    if length(StartPhase)<length(PhaseMask)
        StartPhase(length(PhaseMask))=0;
    end
   ZCF_Best=StartPhase(PhaseMask>0); 
end


%optimize
Options=optimset;
Options.MaxIter=5000;
Options.Display='iter';
[ZCF_Best,Fval] = fminsearch(@costCRLB,double(gather(ZCF_Best)),Options,P,Photons,Bg,PhaseMask);

PFound=P;

Phase=PhaseMask;
Phase(PhaseMask>0)=ZCF_Best;
PFound.ZC_Phase=Phase;

[PSF,PFound]=smi_psf.PointSpreadFunction.scalarPSFZernike(PFound); 
[CRLB]=smi_psf.PointSpreadFunction.crlbPSFPupil(PFound,Photons,Bg)
PSF=gather(PSF);
dipshow(PSF)

PSFStruct=PFound;
end

function out=costCRLB(ZC,PSFStruct,Photons,Bg,PhaseMask)
%to find Pupil:

PlotFlag=0;

%Fill in phase in correct places. 
if length(PSFStruct.ZC_Phase)>length(PhaseMask) %use these
    Phase=PSFStruct.ZC_Phase;
else
    Phase=PhaseMask;
end

Phase(PhaseMask>0)=ZC;

PSFStruct.ZC_Phase=Phase;

[PSF,PSFStruct]=smi_psf.PointSpreadFunction.scalarPSFZernike(PSFStruct);
[CRLB,DET]=smi_psf.PointSpreadFunction.crlbPSFPupil(PSFStruct,Photons,Bg,PlotFlag);

out=double(gather(sum(DET(:))));


end
