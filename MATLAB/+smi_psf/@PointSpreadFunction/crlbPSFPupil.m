function  [CRLB,DET]=crlbPSFPupil(PSFStruct,Photons,Bg,PlotFlag)
%crlbPSFPupil Cramer Rao Lower Bound on (Y,X,Z,Photons,Bg)
%   The PSF is generated using the following model:
%   PSF = |F[OTF]|^2
%   % more...
%
%   PSF is normalized such that the integral over all space = 1
%
% INPUTS:
%   PSFStruct:  PSF Structure.  (Default = createPSFStruct())
%   Photons:    Integrated Emitter Photons. Scalar or Nx1 where N is length(Z)
%   Bg:         Background Photons per Pixel. Scalar or Nx1 where N is length(Z)
%   PlotFlag:   Show plot 0 or 1.  (Default = 1)
%
% OUTPUTS:
%   CRLB        Cramer Rao Lower Bound
%   DET         Determinant of the inverse Fisher Information matrix
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%

%Check for GPU, if not present, overwrite gpuArray
% try
%     gpuDevice
% catch
%     gpuArray=@(x)x;
% end

if nargin <1
    P=smi_psf.PointSpreadFunction.createPSFStruct();
else
    P=PSFStruct;
end

if nargin <2
    Photons=1000;
end

if nargin <3
    Bg=20;
end

if nargin <4
    PlotFlag=1;
end

if length(Photons)==1
   Photons=Photons*ones(length(P.Z),1); 
end

if length(Bg)==1
   Bg=Bg*ones(length(P.Z),1); 
end

%Build Fisher Information Matrix using Numerical Derivatives

%need d model/d Theta


DZ = .05; %micron  step for derivatives

%loop over z

for nn=1:length(P.Z)
    
    PTest=P;
    PTest.Z=P.Z(nn)+[-DZ,0,DZ];
    PTest.SZ=P.SZ+2;
    ModelNorm=smi_psf.PointSpreadFunction.scalarPSFPupil(PTest);
    Model=Photons(nn)*ModelNorm(2:end-1,2:end-1,2)+Bg(nn);
    
    %Derivative in Y. Numerical using whole pixel steps
    D(:,:,1) = Photons(nn)*(ModelNorm(3:end,2:end-1,2)-ModelNorm(1:end-2,2:end-1,2))/(2*P.PixelSize);
    D(:,:,2) = Photons(nn)*(ModelNorm(2:end-1,3:end,2)-ModelNorm(2:end-1,1:end-2,2))/(2*P.PixelSize);
    D(:,:,3) = Photons(nn)*(ModelNorm(2:end-1,2:end-1,3)-ModelNorm(2:end-1,2:end-1,1))/(2*DZ);
    D(:,:,4) = ModelNorm(2:end-1,2:end-1,2);
    D(:,:,5) = Model*0+1;
    %dipshow(gather(D))
    %Build FI
    FI=gpuArray(zeros(5,'single'));
    ModelInv=1./(Model+.1e-5);%dipshow(gather(ModelInv))
    for ii=1:5
        for jj=ii:5
            A=D(:,:,ii).*D(:,:,jj).*ModelInv;
            FI(ii,jj)=sum(A(:));
        end
    end
    
    for ii=1:5
        for jj=ii+1:5
            FI(jj,ii)=FI(ii,jj);
        end
    end
    
    invFI=inv(FI);
    CRLB(nn,:)=diag(invFI);
    DET(nn,:)=det(invFI(1:3,1:3));
end

if PlotFlag
    %figure(1122)
    figure
    plot(P.Z,sqrt(CRLB(:,1))*1000,'R','linewidth',2)
    hold on
    plot(P.Z,sqrt(CRLB(:,2))*1000,'G','linewidth',1.5)
    plot(P.Z,sqrt(CRLB(:,3))*1000,'B','linewidth',1.5)
    xlabel('Z position')
    ylabel('sqrt(crlb) (nm)')
    legend('Y','X','Z')
    axis([min(P.Z) max(P.Z) 0 100])
    hold off
    drawnow;
end

end
