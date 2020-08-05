function [Model,Data] = gaussBlobROIStack(SZ,SMD,VarianceIm,Covariance,PixType)
%gaussBlobROIStack Generates a stack of 2D images containing a single blob
% 
%   A stack of 2D images is generated with a Guassian blob in each image.
%   A blob centered at 0,0 is the center of the upper-left pixel
%
% INPUTS:
%   SZ:         Box size in pixels
%   SMD:        SMD Structre with fields: 
%       PSFSigma:   Gaussian Sigma in Pixels.  Scalar, 1x2, Nx1 or Nx2 [Y X].  
%       Photons:    Integrated Photons in Blob. Scalar or Nx1. 
%       Bg:         Background (Photons per Pixel). Scalar or Nx1. 
%       X:          Center Location of Blob (Pixels)
%       Y:          Center Location of Blob (Pixels)
%   VarianceIm:     Additional Pixel-Dependent ReadNoise Image given as variance (Default=zeros(SZ,SZ))
%   Covariance: Only used when 'Sample' is selected as PixType
%   PixType:    Method of modeling pixel value, can be 'sample' in which
%               case the pixel value is sampled from a gaussian distribution at
%               the center of the pixel or 'integrate' in which case the
%               pixel value is the integrated value of a gaussian
%               distribution over the pixel size. (Default = 'integrate')
%
% OUTPUTS:
%   Model:      Noise Free Image with background as offset
%   Data:       Poisson Noise Corrupted Images with Readnoise (Optional)
%
% REQUIRES:
%   Statistics and Machine Learning Toolbox
%   Parallel Processing Toolbox
%   NVidia GPU
%   cuda_gaussBlobROIStack.ptx
%   cuda_gaussBlobROIStack.cu
% 
% CITATION:
%   Marjolein Meddens 2017, Lidke Lab

%Check input
if nargin<1
    SZ=8;
end
if nargin <2
    SMD=[];
end


if isempty(SMD)
    NFrames=1000;
    SMD.X=SZ/2*ones(NFrames,1);
    SMD.Y=SZ/2*ones(NFrames,1);
    SMD.Photons=1000*ones(NFrames,1);
    SMD.Bg=5*ones(NFrames,1);
    SMD.PSFSigma=1.3*ones(NFrames,1);
else
    NFrames=length(SMD.X);
end

if nargin<3
    VarianceIm=zeros(SZ,SZ,'single');
end
if nargin <4
    Covariance = 0;
end
if nargin <5
    PixType = 'integrate';
end

% expand scalar inputs
if size(SMD.PSFSigma,2) == 1
    SMD.PSFSigma=[SMD.PSFSigma SMD.PSFSigma];
end
if size(SMD.PSFSigma,1) == 1
    SMD.PSFSigma=repmat(SMD.PSFSigma,[NFrames 2]);
end
if size(SMD.Photons,1) == 1
    SMD.Photons=repmat(SMD.Photons,[NFrames 1]);
end
if size(SMD.Bg,1) == 1
    SMD.Bg=repmat(SMD.Bg,[NFrames 1]);
end
if size(Covariance,1) == 1
    Covariance=repmat(Covariance,[NFrames 1]);
end


%%
% allocate output
Model = zeros(SZ,SZ,NFrames,'single');

% create cuda kernels
k_sample = parallel.gpu.CUDAKernel('cuda_gaussBlobROIStack.ptx','cuda_gaussBlobROIStack.cu','kernel_guassiansampleblobs');
k_integ = parallel.gpu.CUDAKernel('cuda_gaussBlobROIStack.ptx','cuda_gaussBlobROIStack.cu','kernel_guassianintegrateblobs');

%gpuDevice gives GPU hardware info
g = gpuDevice;
%Find how many loops need to be run so that ROIimStack fits in GPU memory
Nelem = SZ*SZ*NFrames;
Nloops = ceil(4*4*Nelem/(g.TotalMemory));
%number of frames in each chunk.
NFramesChunk = ceil(NFrames/Nloops);

%The loop sends one chunk of data at each iteration to the gpu.
for ii = 1:Nloops
    if ii == Nloops
        NFramesThisChunk = NFrames - (NFramesChunk*(Nloops-1));
    else
        NFramesThisChunk = NFramesChunk;
    end
    % crop data for this chunk
    startFr = ((ii-1)*NFramesChunk)+1;
    endFr = ii*NFramesChunk;
    endFr = min(endFr,NFrames);
    Ph = single(SMD.Photons(startFr:endFr));
    Bg = single(SMD.Bg(startFr:endFr));
    xCenters = single(SMD.X(startFr:endFr));
    yCenters = single(SMD.Y(startFr:endFr));
    xSigma = single(SMD.PSFSigma(startFr:endFr,2));
    ySigma = single(SMD.PSFSigma(startFr:endFr,1));
    Covar = single(Covariance(startFr:endFr));
    % allocate output variables
    subStack = zeros(SZ,SZ,NFramesThisChunk,'single');
    
    switch PixType
        case 'sample'
            % calculate grid size, fill at least as many blocks as there are processors
            if NFramesThisChunk < g.MultiprocessorCount
                NumBlocks = NFramesThisChunk;
            elseif NFramesThisChunk > (g.MultiprocessorCount * k_sample.MaxThreadsPerBlock)
                NumBlocks = ceil(NFramesThisChunk / k_sample.MaxThreadsPerBlock);
            else
                NumBlocks = g.MultiprocessorCount;
            end
            % calculate block size, multiples of 32 are most efficient
            NumThreads = ceil((NFramesThisChunk/32)/NumBlocks)*32;
            % Setting up number of threads and blocks
            k_sample.GridSize(1) = NumBlocks;
            k_sample.ThreadBlockSize(1) = NumThreads;
            % run cuda function
            [subStack] = feval(k_sample,int32(SZ),int32(NFramesThisChunk),xCenters,yCenters,...
                Ph,Bg,xSigma,ySigma,Covar,   subStack);
        case 'integrate'
            % calculate grid size, fill at least as many blocks as there are processors
            if NFramesThisChunk < g.MultiprocessorCount
                NumBlocks = NFramesThisChunk;
            elseif NFramesThisChunk > (g.MultiprocessorCount * k_integ.MaxThreadsPerBlock)
                NumBlocks = ceil(NFramesThisChunk / k_integ.MaxThreadsPerBlock);
            else
                NumBlocks = g.MultiprocessorCount;
            end
            % calculate block size, multiples of 32 are most efficient
            NumThreads = ceil((NFramesThisChunk/32)/NumBlocks)*32;
            % Setting up number of threads and blocks
            k_integ.GridSize(1) = NumBlocks;
            k_integ.ThreadBlockSize(1) = NumThreads;
            % run cuda function
            [subStack] = feval(k_integ,int32(SZ),int32(NFramesThisChunk),...
                xCenters,yCenters,Ph,Bg,xSigma,ySigma,   subStack);
    end
    % retrieve data from gpu
    Model(:,:,startFr:endFr) = gather(subStack);
end

if nargout > 1
    % add poisson noise
    Data = poissrnd(Model); % this is the slowest step, will need gpu implementation
    % add readnoise
    NoiseIm = sqrt(VarianceIm);
    Data = Data + randn(size(Data)).*repmat(NoiseIm,[1 1 NFrames]);
end

end

