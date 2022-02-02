function [GaussIm] = gaussianImage(SMD,SRImageZoom,ScalebarLength)
%gaussianImage creates gaussian-blob image from SR data.
%
% INPUTS:
%   SMR              SMD results structure with fields X, Y, Bg, X_SE, Y_SE,
%                    Photons, FrameNum, XSize, YSize
%   SRImageZoom      zoom factor [default value of 10]
%   ScalebarLength   scalebar length (um) [default: 10 um]
%
% OUTPUT:
%   GaussIm - Gaussian blob image 
%
% REQUIRES:
%   smi_sim.GaussBlobs.gaussBlobImage.m (and requirements within)
%   Matlab

% Created by:
%    Sandeep Pallikkuth, Lidke Lab 2017

% checking inputs
if nargin<1
    error('Not enough input. Please input SMR structure and SRImageZoom');
elseif nargin<2
    disp('No input on SRImageZoom. Default value of 10 assigned');
    SRImageZoom=10;
end

if ~exist('ScalebarLength', 'var')
   ScalebarLength = 10; % um
end
% creating inputs for smi_sim.GaussBlobs.gaussBlobImage
RawImageSize=[SMD.YSize SMD.XSize];
SZ=RawImageSize*SRImageZoom; %since we need big enough box sizes
SMDin=SMD;
SMDin.XSize=SZ(1);
SMDin.YSize=SZ(2);
SMDin.Bg=zeros(size(SMD.Bg));
SMDin.X=SMD.X*SRImageZoom;
SMDin.Y=SMD.Y*SRImageZoom;
SMDin.PSFSigma=[SMD.Y_SE,SMD.X_SE]*SRImageZoom;
SMDin.Photons=ones(size(SMD.Photons));
SMDin.FrameNum=ones(size(SMD.FrameNum));
SMDin.NFrames=1;
% creating GaussIm
GaussIm = smi_sim.GaussBlobs.gaussBlobImage(SMDin);
GaussIm = smi.BaGoL.scaleIm(GaussIm);
GaussIm = smi_vis.GenerateImages.scalebar(GaussIm, SMD.PixelSize, ...
                                          ScalebarLength);
GaussIm = smi_vis.GenerateImages.colorImage(GaussIm);

end
