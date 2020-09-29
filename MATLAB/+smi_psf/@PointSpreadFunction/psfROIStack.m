function [Model,Data]=psfROIStack(PSF,XYSamPerPix,ZSamPerUnit,SZ,SMD,NoiseIm) 
%psfROIStack generates boxes of the data containing single particles (blobs).
%
%psfROIStack() gets the size of the boxes, the structure SMD, which
%contains the X, Y, Z positions, photons and the background of the boxes,
%PSF, which is the sampled PSF,XYSamPerPix and ZSamPerUnit, which are the 
%number of the samples in PSF along X,Y ans Z axes, and finally NoiseIm is
%the read-out noise and generates boxes of the data containing single
%particles (blobs).
%
% INPUTS
%   PSF:        The sampled PSF (Point Spread Function).
%   XYSamPerPix:The number of the PSF samples along the X and Y axes per
%               pixel (No default).       
%   ZSamPerUnit:The number of the samples along the Z-axis. (No default)
%   SZ:         size of the box (Default=15) (Pixels).
%   SMD:        This is a structure containing some fileds. We use the fields 
%               X Y, Z, Photons and Background.
%   NoiseIm:    The read-out noise. (Defualt = 0), (Photons).
%
% OUTPUTS:
%   Model:      The box with the blob, which is not croppted with the
%               shot-noise.
%   Data:       The model corrupted with the shot-noise (Poisson noise).
%
% REQUIRES:
%   Parallel processing toolbox
%   NVidia GPU
%   psfSample3DBlob.cu
%   psfSample3DBlob.ptx

% Created by:
%   Mohamadreza Fazel (Lidke Lab, 2017)
%

if nargin < 3
   error('The Sampled PSF and the number of samples along the different axes must be given.'); 
end
if nargin < 4
   SZ = 15;
end
if nargin < 5
   SMD.X = rand([1,10000])*SZ/3+SZ/3;
   SMD.Y = rand([1,10000])*SZ/3+SZ/3;
   SMD.Z = (rand([1,10000])*size(PSF,3)/ZSamPerPix)-size(PSF,3)/(ZSamPerPix*2);
   SMD.Photons = 1000*ones([1,10000]);
   SMD.Bg = 15;
end
if nargin < 6
   NoiseIm = zeros(SZ); 
end
[a,b,c]=size(PSF);
Filter = ones(XYSamPerPix);
PSFIntegral = zeros(a+XYSamPerPix-1,b+XYSamPerPix-1,c);
for ii=1:c
    PSFIntegral(:,:,ii) = (1/XYSamPerPix)^2*conv2(PSF(:,:,ii),Filter);
end
%After the convolution, some of the resulted componenets at the begining and 
%at the end of the array do not represent a whole pixel so we delete them.
PSFIntegral = PSFIntegral(XYSamPerPix:end-(XYSamPerPix-1),XYSamPerPix:end-(XYSamPerPix-1),:);

X = SMD.X;
Y = SMD.Y;
Z = SMD.Z;
Photons = SMD.Photons;
if isscalar(SMD.Bg)
    Bg = SMD.Bg*ones(size(SMD.X));
else
    Bg = SMD.Bg;
end
FrameNum = length(X);
Model = zeros(SZ,SZ,FrameNum,'single');
Nelem =SZ*SZ*FrameNum;
kc = parallel.gpu.CUDAKernel('cuda_PSFSample3DBlob.ptx','cuda_PSFSample3DBlob.cu','cuda_PSFSample3DBlob');
%gpuDevice gives the info of your hardware, like the available
%memory (TotalMemory).
g = gpuDevice;
%The whole set of data cannot be sent into the gpu at the same
%time because of the limitation of the gpu memory. Here we are
%finding the number of chuncks of data based on the available
%memory (g.Totalmemory) to send them to the gpu.
Nloops = ceil(4*4*Nelem/(g.TotalMemory));
%number of frames in each chunk.
NFrameChunk = floor(FrameNum/Nloops);
Cent = round(size(PSFIntegral,1)/2);

for ii = 1:Nloops
       
       if ii == Nloops
          NFrameThisChunk = FrameNum - (NFrameChunk*(Nloops-1));
       else
          NFrameThisChunk = NFrameChunk;
       end
       EndFrame = (ii-1)*NFrameChunk+NFrameThisChunk;
       StartFrame = (ii - 1)*NFrameThisChunk + 1;
       %The chunk of data that are sent to the gpu
       ThisData = zeros(SZ,SZ,NFrameThisChunk);
       %computing the number of the blocks.
       %NumBlocks = ceil(NFrameThisChunk / kc.MaxThreadsPerBlock);
       if NFrameThisChunk < g.MultiprocessorCount
           NumBlocks = NFrameThisChunk;
       elseif NFrameThisChunk > (g.MultiprocessorCount * kc.MaxThreadsPerBlock)
           NumBlocks = ceil(NFrameThisChunk / kc.MaxThreadsPerBlock);
       else
           NumBlocks = g.MultiprocessorCount;
       end
       %calculating the number of the threads in each block. It is more
       %efficient when it is proportional to 32.
       NumThreads = ceil(NFrameThisChunk/(32*NumBlocks))*32;
       %setting the number of the threads in each block
       kc.ThreadBlockSize(1) = NumThreads;
       %setting the number of the blocks
       kc.GridSize(1) = NumBlocks;
       Bg_in = zeros(1,NFrameThisChunk,'single');
       %calling cuda-code

       [ThisData] = feval(kc,X,Y,Z,Photons,Bg_in,XYSamPerPix,ZSamPerUnit,int32(SZ),int32(Cent),int32(size(PSFIntegral)),...
           PSFIntegral,NFrameThisChunk,ThisData);
       %retrieve output from gpu memory
       Model(:,:,StartFrame:EndFrame) = gather(ThisData);
end
ModelCoef = sum(sum(Model));
ModelCoef = ModelCoef(:);
for nn = 1:size(Model,3)
    Model(:,:,nn) = Photons(nn)*Model(:,:,nn)/ModelCoef(nn) + Bg(nn);
end
if nargout > 1
    % add poisson noise
     Data = poissrnd(Model); 
    % add readnoise
    Data = Data + randn(size(Data)).*repmat(NoiseIm,[1 1 FrameNum]);
end
end
