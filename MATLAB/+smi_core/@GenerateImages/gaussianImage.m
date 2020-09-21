function [GaussIm] = gaussianImage(SMR,SRImageZoom)
% Function gaussianImage creates gaussian-blob image from SR data
% INPUT:
%   SMR - SMR results structure with fields X, Y, Bg, X_SE, Y_SE, Photons,
%   FrameNum, XSize, YSize
%   SRImageZoom - zoom factor (default value of 10)
% OUTPUT:
%   GaussIm - Gaussian blob image 
%
% REQUIRES:
%   smi_sim.GaussBlobs.gaussBlobImage.m (and requirements within)
%   Matlab
%
% CITATION:
%   Sandeep Pallikkuth, Lidke Lab 2017

% checking inputs
if nargin<1
    error('Not enough input. Please input SMR structure and SRImageZoom');
elseif nargin<2
    disp('No input on SRImageZoom. Default value of 10 assigned');
    SRImageZoom=10;
end

% creating inputs for smi_sim.GaussBlobs.gaussBlobImage
RawImageSize=[SMR.YSize SMR.XSize];
SZ=RawImageSize*SRImageZoom; %since we need big enough box sizes
SMDin.Bg=zeros(size(SMR.Bg));
SMDin.X=SMR.X*SRImageZoom;
SMDin.Y=SMR.Y*SRImageZoom;
SMDin.PSFSigma=[SMR.Y_SE,SMR.X_SE]*SRImageZoom;
SMDin.Photons=ones(size(SMR.Photons))*1000;
SMDin.FrameNum=ones(size(SMR.FrameNum));
% creating GaussIm
GaussIm = smi_sim.GaussBlobs.gaussBlobImage(SZ,1,SMDin);

end
