function [PSF,PSFStruct_OS] = oversamplePFSPupil(PSFStruct,Sampling)
%oversamplePFSPupil Resamples a PSF to smaller pixels.  
%   Detailed explanation goes here. Yes, detailed. All help 
%   must have the same format. This tempate describes what is required in
%   the help and also defines a formatting convention.  Check the
%   formatting of your function using 'doc SMA_Core', 
%   'help SMA_Core.helpTemplate' and 'doc SMA_Core.helpTemplate'
%  
% INPUTS:
%   PSFStruct:      PSF Structure. Must contain 'Pupil' field. 
%   Sampling:       Oversampling. Bust be integer >=1 (Default=4)
%
% OUTPUTS:
%   PSF:            Resampled PSF
%   PSFStruct_OS:   Updated PSF Structure
%
% REQUIRES:
%   Statistics Toolbox
%   Parallel Procesing Toolbox
%   NVidia GPU
%
% CITATION:
%   Full citatation with all authors.         


if nargin<2
    Sampling=4;
end

P_OS=PSFStruct;
P_OS.SZ=PSFStruct.SZ*Sampling;
P_OS.PixelSize=PSFStruct.PixelSize/Sampling;

%Pad pupil
Pupil_1X=gather(P_OS.Pupil);
PMag=extend(Pupil_1X(:,:,1),[Sampling,Sampling]*P_OS.OSZ,'symmetric');
PPhase=extend(Pupil_1X(:,:,2),[Sampling,Sampling]*P_OS.OSZ,'symmetric');
P_OS.Pupil=single(cat(3,PMag,PPhase));

%update OTF size
P_OS.OSZ=Sampling*PSFStruct.OSZ;

PSF=gather(smi_psf.PointSpreadFunction.scalarPSFPupil(P_OS));
PSFStruct_OS=P_OS;

end
