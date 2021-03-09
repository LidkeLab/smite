classdef FindROI < handle
    % FindROI Finds and collates subregions from a stack of 2D images
    %
    % Subregions are found by filtering with a difference of Gaussians filter
    % followed by finding the local maximum.  Local maxima are used as the
    % center of ROIs if the estimated single molecule intesity is greater than
    % MinPhotons.
    %
    % Gaussian filtering is done using the recursive method of Young.
    % Local maximum finding is done using a novel separable filter (Lidke).
    %
    % Filtering and maximum finding are implemented on GPU using CUDA
    % and the Parallel Processsing Toolbox.
    %
    %
    % REQUIRES:
    %   MATLAB 2014a or later versions
    %   Parallel Procesing Toolbox
    %   NVidia GPU
    %   smi_cuda_FindROI.ptx
    %   smi_cuda_FindROI.cu
    %
    % CITATION:
    %   Ian T. Young, Lucas J. van Vliet,
    %   Recursive implementation of the Gaussian filter,
    %   Signal Processing, Volume 44, Issue 2, 1995
    
    
    properties
        Data            %Stack of 2D images
        BoxSize=7       %Linear box size for fitting (Pixels)(Default=7)
        BoxOverlap=2    %Overlap of boxes (Integer Pixels)(Default=2)
        MinPhotons=200  %Minimum number of photons from emitter (Default=200)
        PSFSigma=1.3    %Sigma of 2D Gaussian PSF (Pixels) (Default=1.3)
        ROIs            %Stack of subregions
        LocalMaxIm      %Binary Image showing local maxima above the threshold
        PlotBoxFrame=1  %If Verbose >= 3, plot boxes for this frame (Default=1)
        Verbose=1       %Verbosity level
    end
    
    properties (Access = protected)
        SigmaSmall      %Small smoothing kernel (Pixels)(Default=PSFSigma)
        SigmaLarge      %Large smoothing kernel (Pixels)(Default=2*PSFSigma) 
    end
    
    methods
        function [obj]=FindROI(SMF,Data)
        %FindROI Create a FindROI object
        %[obj]=FindROI(SMF,Data)
        %
        % [obj]=FindROI()
        % Creates an object with default property values.
        %
        % [obj]=FindROI(SMF)
        % Creates an object and copies SMF fields to Properties
        %
        % [obj]=FindROI(SMF,Data)
        % Creates an object and copies Data and SMF fields to Properties
        %
        % INPUTS
        %   SMF:    SMF structure with fields (Optional)
        %       BoxFinding:
        %           BoxSize         %Linear box size for fitting (Pixels)
        %           BoxOverlap      %Overlap of boxes allowed (Pixels)
        %           MinPhotons      %Minimum number of photons from emitter
        %       Fitting:
        %           PSFSigma        %Sigma of 2D Gaussian PSF Model (Pixels)
        %   Data:   Stack of 2D images (Optional)
        %
        % OUTPUTS
        %   obj:    FindROI object
        %
            if nargin > 0
                obj.BoxSize=SMF.BoxFinding.BoxSize;
                obj.BoxOverlap=SMF.BoxFinding.BoxOverlap;
                obj.MinPhotons=SMF.BoxFinding.MinPhotons;
                obj.PSFSigma=mean(SMF.Fitting.PSFSigma);
                obj.SigmaSmall=obj.PSFSigma;
                obj.SigmaLarge=2*obj.PSFSigma;
            end
            
            if nargin > 1
                obj.Data=Data;
            end

        end
        
        function [ROIStack,SMD]=findROI(obj,SMD)
            %findROI Main method of FindROI
            %[ROIs]=obj.findROI()
            %
            % Finds and cuts out subregions of size BoxSize x BoxSize. If no
            % SMD is given as input, a default SMD is created.  The XBoxCorner
            % and YBoxCorner are modified in the SMD structure.
            %
            % INPUTS
            %   SMD      SMD data structure (Optional)
            %
            % OUTPUTS
            %   ROIStack Stack of N 2D subregions (BoxSize x BoxSize x N)
            %   SMD      SMD data structure based on input SMD structure if
            %            provided, or newly created.  Fields defined here are:
            %               XBoxCorner, YBoxCorner, FrameNum, NFrames, XSize,
            %               YSize, NDims (findROI assumes 2D)
            
            %Create SMD if needed
            if nargin<2
               SMD=smi_core.SingleMoleculeData.createSMD();
            end
            
            %Calculate min allowed local maximum value
            SigmaBlob=sqrt(obj.SigmaSmall^2+mean(obj.PSFSigma)^2);
            MinVal = obj.MinPhotons*(1/(2*pi*SigmaBlob^2));
            
            %Calculate Kernel Size
            LMKernelSize=obj.BoxSize-obj.BoxOverlap;
            
            %Break data into chunks that fit in GPU memory
            NCopies=6; %for out of place operations
            NBytesPerPixel=4; %single float
            NElem = numel(obj.Data);
            g = gpuDevice;
            NLoops = ceil(NBytesPerPixel*NCopies*NElem/(g.AvailableMemory));
            
            %Number of frames in each chunk.
            ZSize = size(obj.Data,3);
            ZChunkSize = floor(ZSize/NLoops);
            LocalMaxIm_Temp = zeros(size(obj.Data));
            
            %The loop sends one chunk of data at each iteration to the gpu.
            for ii = 1:NLoops
                
                ZStart=(ii-1)*ZChunkSize+1;
                ZEnd = ii*ZChunkSize;
                if ii == NLoops
                    ZEnd = ZSize;
                end
                
                %The chunk of data that is sent to the GPU
                SubStack = obj.Data(:,:,ZStart:ZEnd);
                
                %Call static methods
                [D_A] = smi_core.FindROI.gaussInPlace(SubStack,obj.SigmaSmall);
                [D_B] = smi_core.FindROI.gaussInPlace(SubStack,obj.SigmaLarge);
                D_A=D_A-D_B;
                [SubLocalMax]=smi_core.FindROI.localMax(D_A,LMKernelSize,MinVal);
                
                %Collect Results
                LocalMaxIm_Temp(:,:,ZStart:ZEnd) = SubLocalMax;
            end
            wait(g);
            obj.LocalMaxIm=LocalMaxIm_Temp;
            
            %Finding coordinates of the local maxima
            XSize = size(obj.Data,2);
            YSize = size(obj.Data,1);
            
            [Ind] = find(obj.LocalMaxIm);
            [YInd,XInd,ZInd]=ind2sub([YSize XSize ZSize],Ind);
            
            %Finding the boxes
            HalfBox = floor(obj.BoxSize/2+0.5);
            StartRow = YInd-HalfBox+1;
            minusRow = find(StartRow<1);
            StartRow(minusRow) = 1;
            StartCol = XInd-HalfBox+1;
            MinusCol = find(StartCol<1);
            StartCol(MinusCol) = 1;
            
            EndRow = StartRow + obj.BoxSize-1;
            EndRow(minusRow) = obj.BoxSize;
            LargeRow = find(EndRow > YSize);
            EndRow(LargeRow) = YSize;
            StartRow(LargeRow) = YSize - obj.BoxSize +1;
            
            EndCol = StartCol + obj.BoxSize-1;
            EndCol(MinusCol) = obj.BoxSize;
            LargeCol = find(EndCol > XSize);
            EndCol(LargeCol) = XSize;
            StartCol(LargeCol) = XSize - obj.BoxSize + 1;
            
            %Outputs:
            N = length(StartRow); 

            ROIStack = zeros(obj.BoxSize,obj.BoxSize,N);
            for nn = 1:N
                ROIStack(:,:,nn) = obj.Data(StartRow(nn):EndRow(nn), StartCol(nn):EndCol(nn), ZInd(nn));
            end
            SMD.XBoxCorner = StartCol;
            SMD.YBoxCorner = StartRow;
            SMD.FrameNum = ZInd;
            SMD.NFrames = ZSize;
            SMD.XSize = XSize;
            SMD.YSize = YSize;
            SMD.NDims = 2;

            if Verbose >= 3
                plotBox(SMD, obj.Data, obj.PlotBoxFrame, obj.BoxSize)
            end
        end
        
        function showOverlay(obj)
        %showOverlay Show local maxima on top of Data
        % 
        % REQUIRES:
        %   DipImage
        
        overlay(stretch(obj.Data),obj.LocalMaxIm)
            
        end
        
    end
    
    methods (Static)
        function [Stack] = gaussInPlace(Stack, Sigma)
        %gaussInPlace() In-place Gaussian filter  
        % [SubStack] = gaussInPlace(SubStack, Sigma)
        %
        % Filtering is applied using the recursive method of Young and Van
        % Vliet and implimented on the GPU using CUDA. 
        %
        % This is an in-place filter, meaning that input array is
        % modified and given back as the ouput.  
        %
        % If the input is not gpuArray, it is copied to the GPU, which has
        % the effect that the input is not modified. 
        % 
        % INPUTS:
        %   Stack:      Stack of 2D images to be filtered
        %   Sigma:      Gaussian Kernel (Pixels)
        %
        % OUTPUTS:
        %   Stack:      Modified 2D stack of images (gpuArray) 
        %
        % REQUIRES:
        %   MATLAB 2014a or later versions
        %   Parallel Procesing Toolbox
        %   NVidia GPU
        %   smi_cuda_FindROI.ptx
        %   smi_cuda_FindROI.cu
        %
        % CITATION:
        %   Ian T. Young, Lucas J. van Vliet,
        %   Recursive implementation of the Gaussian filter,
        %   Signal Processing, Volume 44, Issue 2, 1995
        
        %Setting up constants for the recursive method (See citation)
        if Sigma > 2.5
            Q = 0.98711*Sigma-0.96330;
        else
            Q = 3.97156-4.14554*sqrt(1-0.26891*Sigma);
        end
        
        Qq = Q*Q;
        Qqq = Qq*Q;
        
        B0 = 1.57825+2.44413*Q+1.4281*Qq+0.422205*Qqq;
        B1 = 2.44413*Q+2.85619*Qq+1.26661*Qqq;
        B2 = -1.4281*Qq-1.26661*Qqq;
        B3 = 0.422205*Qqq;
        B = 1 - (B1+B2+B3)/B0;
        
        %Creating GPU CUDA kernel objects from PTX and CU code
        K_Gauss = parallel.gpu.CUDAKernel('smi_cuda_FindROI.ptx','smi_cuda_FindROI.cu','kernel_gaussMajor');
        
        %Calling the gpu code to apply Gaussian filter along major
        %This is an in-place operation
        K_Gauss.GridSize(1) = size(Stack,3);
        K_Gauss.ThreadBlockSize(1) = size(Stack,2);
        Stack = feval(K_Gauss,Stack,size(Stack,1),B0,B1,B2,B3,B);
        
        %Permuting and doing the other dimension
        Stack=permute(Stack,[2 1 3]);
        
        K_Gauss.GridSize(1) = size(Stack,3);
        K_Gauss.ThreadBlockSize(1) = size(Stack,2);
        Stack = feval(K_Gauss,Stack,size(Stack,1),B0,B1,B2,B3,B);
        
        %Permute back.  
        Stack=permute(Stack,[2 1 3]);
        end
        
        function [LocalMaxIm]=localMax(Stack,KernelSize,MinVal)
        %localMax() Generates a binary image of local maxima
        %
        % Local maxima are found within a square area of size KernelSize+1.
        % This is a 2D operation applied to each image in the Stack.  The
        % output is a binary image that has true pixels at the locations of
        % the local maxima that are larger MinVal. 
        % 
        % Uses a x,y separable, out-of-place filter to find local maxima
        % in 2D images. 
        % 
        % INPUTS: 
        %   Stack:      Stack of 2D images. 
        %   KernelSize: Closest allowed local maxima in each dimension
        %   MinVal:     Minimum value of local maxima to be valid
        %
        % OUTPUTS: 
        %   LocalMaxIm: Stack of binary images indicating local maximum     
        %
        % REQUIRES:
        %   MATLAB 2014a or later versions
        %   Parallel Procesing Toolbox
        %   NVidia GPU
        %   smi_cuda_FindROI.ptx
        %   smi_cuda_FindROI.cu
        
            K_LM1 = parallel.gpu.CUDAKernel('smi_cuda_FindROI.ptx','smi_cuda_FindROI.cu','kernel_LocalMaxFirstPass');
            K_LM2 = parallel.gpu.CUDAKernel('smi_cuda_FindROI.ptx','smi_cuda_FindROI.cu','kernel_LocalMaxSecondPass');
            
            %Make out-of-place storage and final image
            TempStack=gpuArray(Stack);
            LocalMaxIm=gpuArray(Stack);
            
            K_LM1.GridSize = [size(Stack,2),size(Stack,3),1];
            K_LM1.ThreadBlockSize(1) = size(Stack,1);
            [TempStack] = feval(K_LM1,Stack,TempStack,KernelSize,MinVal);
            
            %Permute for faster cuda kernel
            TempStack=permute(TempStack,[2 1 3]);
            LocalMaxIm=permute(LocalMaxIm,[2 1 3]);
             
            K_LM2.GridSize = [size(TempStack,2),size(TempStack,3),1];
            K_LM2.ThreadBlockSize(1) = size(TempStack,1);
            [LocalMaxIm] = feval(K_LM2,TempStack,LocalMaxIm,KernelSize,MinVal);
            
            %Permute back and give host output
            LocalMaxIm=permute(LocalMaxIm,[2 1 3]);
            LocalMaxIm = gather(LocalMaxIm);
        end

        plotBox(SMD, Data, Frame, BoxSize)
    end
end


