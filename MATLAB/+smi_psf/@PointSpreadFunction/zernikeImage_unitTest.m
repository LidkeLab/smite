function [Report]=zernikeImage_unitTest()
%zernikeImage_unitTest Tests zernikeImage functionality.

Report = 0;

%% Check orthogonality

% These should all be zero
NMax=9;
Radius=64;
SZ=256;
for nn=1:NMax
    for kk=nn+1:NMax
        Im1=gather(smi_psf.PointSpreadFunction.zernikeImage(nn,SZ,Radius));
        Im2=gather(smi_psf.PointSpreadFunction.zernikeImage(kk,SZ,Radius));
        IP=sum(sum(Im1.*Im2))/(pi*Radius^2);
        if abs(IP)>1e-2
            fprintf('N1=%d, N2=%d. Inner product: %g \n',nn,kk,IP)
        end
    end
end


%% Check Normalization

% These should all be Em*pi(2N-2);
NMax=9;
Radius=64;
SZ=256;
for nn=1:NMax
    Im1=gather(smi_psf.PointSpreadFunction.zernikeImage(nn,SZ,Radius));
    [N,M]=smi_psf.Zernike.zNoll2NM(nn);
    IP=sum(sum(Im1.*Im1))/(pi*Radius^2);
    if (N==0)&&(M==0);
        Norm=1;
    else
        if M==0
            Em=2;
        else
            Em=1;
        end
        Norm=Em/(2*N+2);
    end
    fprintf('N=%d, Inner product: %g Expected: %g\n',nn,IP,Norm)
end

Report = 1;
