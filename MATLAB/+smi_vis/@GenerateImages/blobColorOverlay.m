function [OverlayImage] = blobColorOverlay(Sequence,SMD,AutoScale)
%blobColorOverlay creates overlay of fitted emitters onto data.
%   
%   Creates color overlay of fitted emitters (green) onto data (red). 
%   smi_sim.GaussBlobs.gaussBlobImage is used to generate the blobs of the
%   fitted emitters. 
%   
% INPUTS:
%   Sequence:   Data stack onto which fitted emitters will be overlayed
%   SMD:        Structure containing at least the following fields:
%       X:          Gaussian Blob Center X (Pixels) 
%       Y:          Gaussian Blob Center Y (Pixels) 
%       PSFSigma:   Gaussian Sigma (Pixels).  Can be scalar or 1x2. 
%       Photons:    Integrated Photons in Blob 
%       Bg:         Photons per Pixel
%       FrameNum:   Frame number of Gaussian Blob Location
%   AutoScale   scale the intensity of the Model based on max intensity of
%               rawdata (0:not change intensity of the model or 1:change it)
%
% OUTPUTS:
%   OverlayImage:   (Optional) Image of fitted emitters (green) overlayed 
%                   onto data (red). If no output is given, overlay image
%                   will be displayed
%
% REQUIRES:
%   Statistics Toolbox
%   Parallel Procesing Toolbox
%   Dipimage toolbox
%   NVidia GPU

% Created by:
%   Marjolein Meddens and Hanieh Mazloom-Farsibaf 2017, Lidke Lab        

if nargin<3
    AutoScale=0;
end 

% generate blob stack for fits
SZ = [size(Sequence,1),size(Sequence,2)];
NFrames = size(Sequence,3);
[Model] = smi_sim.GaussBlobs.gaussBlobImage(SZ,NFrames,SMD);

% scale the intensity of the Model based on the Max(data)
if AutoScale && max(Sequence(:))>max(Model(:))
    Model2=Model*(2*(max(Sequence(:))/max(Model(:))));
else
    Model2=Model;
end
% overlay onto data
OverlayImage=smi_vis.GenerateImages.rgbImage(Sequence,Model2,Model2*0);

if nargout == 0
    smi_vis.GenerateImages.showim(OverlayImage);
end

end

