classdef GaussMLE < handle
    % GaussMLE Maximum Likelihood Estimate of 2D Gaussian blob parameters
    %
    % GaussMLE fits a 2D Gaussian blob model to the each image in a stack of 2D
    % images.  Several fit types are available. The MLE is maximized with an
    % iterative Newton-Raphson method and implemented on GPU using NVIDIA
    % CUDA as described by Smith et al.
    %
    % The input is assumed to be gain and offset corrected so that
    % pixel values indicate effective photon counts. Noise model is 'Poisson',
    % which is suitable for EMCCD cameras used at high gain, or 'SCMOS' which
    % includes a pixel-wise readnoise as described by Huang et al.
    %
    % REQUIRES:
    %   MATLAB 2014a or later versions
    %   Parallel Procesing Toolbox
    %   NVidia GPU
    %   cuda_FindROI.ptx
    %   cuda_FindROI.cu
    %
    % CITATION:
    %   Smith, C., Joseph, N., Rieger, B. et al.
    %   Fast, single-molecule localization that achieves theoretically minimum
    %   uncertainty. Nat Methods 7, 373–375 (2010).
    %   https://doi.org/10.1038/nmeth.1449
    %
    %   Huang, F., Hartwich, T., Rivera-Molina, F. et al. Video-rate nanoscopy
    %   using sCMOS camera–specific single-molecule localization algorithms.
    %   Nat Methods 10, 653–658 (2013).
    %   https://doi.org/10.1038/nmeth.2488
    
    
    properties
        Data                %Stack of 2D subregions (NxNxM)
        FitType='XYNB'      %One of {'XYNB','XYNBS','XYNBSXSY','XYZNB'} (Default='XYNB')
        Iterations=20;      %Newton-Raphson iterations
        CameraType='EMCCD'  %One of {'EMCCD','SCMOS'} (Default='EMCCD')
        PSFSigma=1.3        %Known or initial PSF Sigma (Pixels) Scalar or [SigmaX SigmaY] (Default=1.3)
        VarianceIm      %The variance image used for SCMOS noise model (Optional)
        XBoxCorner      %Location of fit box within VarianceIm (Pixel) (Optional)
        YBoxCorner      %Location of fit box within VarianceIm (Pixel) (Optional)
        ZFitStruct      %Structure for astigmatic fitting with fields:
            Ax      %
            Ay
            Bx
            By
            Gamma
            D
        Statistics      %Structure with compuational performance info
        SMD
    end
    
    methods
        
        function [obj]=GaussMLE(SMF,Data)
        %   %GaussMLE Create a GaussMLE object
        %[obj]=GaussMLE(Data,SMF)
        %
        % [obj]=GaussMLE()
        % Creates an object with default property values.
        %
        % [obj]=GaussMLE(SMF)
        % Creates an object and copies SMF fields to Properties
        %
        % [obj]=GaussMLE(SMF,Data)
        % Creates an object and copies Data and SMF fields to Properties
        %
        % INPUTS
        %   
        %   SMF:    SMF structure with fields (Optional)
        %       Fitting:
        %           PSFSigma: Sigma of 2D Gaussian PSF Model (Pixels)
        %   Data:   Stack of 2D images (Optional)
        %
        % OUTPUTS
        %   obj:    FindROI object
        %
        
            if nargin > 0
                obj.FitType=SMF.Fitting.FitType;
                obj.PSFSigma=SMF.Fitting.PSFSigma;
                obj.Iterations=SMF.Fitting.Iterations;
                obj.VarianceIm=SMF.Data.CameraReadNoise.^2;
                obj.ZFitStruct=SMF.Fitting.ZFitStruct;
            end
            
            if nargin > 1
                obj.Data=Data;
            end
        
        end
        
        
        function [Results, Statistics]=gaussMLE(obj,SMD)
        %gaussMLE Fits 2D stack of images to Gaussian blob model     
        %
        % INPUTS:
        %   SMD:        SMD structure to be modified (Optional)
        % OUTPUTS:
        %   Results     SMD stucture with modified fields:
        %       X:          Found position in X (Mx1)(Pixels)
        %       Y:          Found position in Y (Mx1)(Pixels)
        %       Z:          Found position in Z (Mx1)(microns) ('XYZNB')
        %       Photons:    Found number of photons from emitter (Mx1)     
        %       Bg:         Found background (Photons/Pixel) (Mx1)
        %       PSFSigma    Found PSFSigma (Mx1)or(Mx2)(Pixels) ('XYNBS','XYNBSXSY')
        %       X_SE:           Uncertainly as standard error (Mx1)
        %       Y_SE:           Uncertainly as standard error (Mx1)
        %       Z_SE:           Uncertainly as standard error (Mx1)
        %       Photons_SE:     Uncertainly as standard error (Mx1)
        %       Bg_SE:          Uncertainly as standard error (Mx1)
        %       PSFSigma_SE     Uncertainly as standard error (Mx1)or(Mx2)
        %       LogLikelihood:  Log likelihood at MLE (Mx1)
        %       PValue:         PValue (Mx1)
        %
        % REQUIRES:
        %       Parallel Computing Toolbox
        %       NVidia GPU
        %       cuda_gaussMLEv2.cu
        %       cuda_gaussMLEv2.ptx
        %
        
        if nargin<2
            SMD=smi_core.SingleMoleculeData.createSMD();
        end
        Results=SMD;
        
        %Set up kernels and number of parameters
        switch obj.FitType
            case 'XYNB'
                NP=4;
                KernelID='_XYNB_';
            case 'XYNBS'
                NP=5;
                KernelID='_XYNBS_';
            case 'XYZNB'
                NP=5;
                KernelID='_XYNBZ_';
            case 'XYNBSXSY'
                NP=6;
                KernelID='_XYNBSXSY_';
        end
        
        switch obj.CameraType
            case 'SCMOS'
                KernelID=['_SCMOS' KernelID(2:end)];
        end
        
        k = parallel.gpu.CUDAKernel('cuda_gaussMLEv2.ptx','cuda_gaussMLEv2.cu',KernelID);
        
        N=size(obj.Data,3);
        SZ=size(obj.Data,1);
        G=gpuDevice();
        AM=G.AvailableMemory;
        BytesPerFloat=4;
        MemoryPerFit=BytesPerFloat*(SZ*SZ+NP+NP+1); %Data plus outputs
        MaxFits=floor(AM/MemoryPerFit/2); %Factor of 2 to be conservative
        
        BSZ=128; %Threads per block hard-coded into kernel
        NKernelsCalls = ceil(N/MaxFits);
        NFitsPerCall=min(MaxFits,N);
        
        %Setup output arrays
        Params_out=zeros(N,NP,'single');
        CRLB_out=zeros(N,NP,'single');
        LL_out=zeros(N,1,'single');
        
        tic
        for nn=1:NKernelsCalls
            StartIndex=(nn-1)*NFitsPerCall+1;
            EndIndex=min((nn+1)*NFitsPerCall,N); %Don't go past data size
            NFitsActual = EndIndex-StartIndex+1;
            SubData=obj.Data(:,:,StartIndex:EndIndex);
            k.GridSize = [ceil(NFitsActual/BSZ) 1];
            k.ThreadBlockSize = [BSZ 1];
            d_Parameters=zeros(NFitsActual,NP,'single');
            d_CRLBs=zeros(NFitsActual,NP,'single');
            d_LogLikelihood=zeros(NFitsActual,1,'single');
            switch obj.CameraType
                case 'EMCCD'
                    switch obj.FitType
                        case {'XYNB','XYNBS','XYNBSXSY'}
                            [P, CRLB,LL] = feval(k,SubData,mean(obj.PSFSigma),SZ,obj.Iterations,d_Parameters,d_CRLBs,d_LogLikelihood,NFitsActual);
                        case 'XYZNB'
                            [P, CRLB,LL] = feval(k,SubData,obj.PSFSigma(1),...
                                obj.ZFitStruct.Ax,obj.ZFitStruct.Ay, obj.ZFitStruct.Bx, obj.ZFitStruct.By,obj.ZFitStruct.Gamma,obj.ZFitStruct.D,obj.PSFSigma,...
                                SZ,obj.Iterations,d_Parameters,d_CRLBs,d_LogLikelihood,NFitsActual);
                        otherwise
                            error('gaussMLE: Unknown fit type: %s',obj.FitType)
                    end
                    
                case 'SCMOS'
                    switch obj.FitType
                        case {'XYNB','XYNBS','XYNBSXSY'}
                            [P, CRLB,LL] = feval(k,SubData,BoxCorners,VarianceImage,...
                                mean(Sigma),SZ,size(VarianceImage,1),obj.Iterations,...
                                d_Parameters,d_CRLBs,d_LogLikelihood,NFitsActual);
                        case 'XYZNB'
                            Z0=zeros(NFitsActual,1,'single');
                            [P, CRLB,LL] = feval(k,SubData,BoxCorners,obj.VarianceIm,Z0,...
                                obj.PSFSigma(1),obj.ZFitStruct.Ax,obj.ZFitStruct.Ay, obj.ZFitStruct.Bx,obj.ZFitStruct.By,...
                                obj.ZFitStruct.Gamma,obj.ZFitStruct.D,obj.PSFSigma(2),...
                                SZ,size(VarianceImage,1),obj.Iterations,d_Parameters,d_CRLBs,d_LogLikelihood,NFitsActual);
                        otherwise
                            error('gaussMLE: Unknown fit type: %s',obj.FitType)
                    end
                otherwise
                    error('gaussMLE: Unknown camera type: %s',obj.CameraType)
            end
            Params_out(StartIndex:EndIndex,:)=gather(P);
            CRLB_out(StartIndex:EndIndex,:)=gather(CRLB);
            LL_out(StartIndex:EndIndex)=gather(LL);
        end
        
        Results.Y=Params_out(:,1);
        Results.X=Params_out(:,2);
        Results.Photons=Params_out(:,3);
        Results.Bg=Params_out(:,4);
        Results.Y_SE=sqrt(CRLB_out(:,1));
        Results.X_SE=sqrt(CRLB_out(:,2));
        Results.Photons_SE=sqrt(CRLB_out(:,3));
        Results.Bg_SE=sqrt(CRLB_out(:,4));
        
        switch obj.FitType
            case 'XYNB'
            case 'XYNBS'
                Results.Sigma=Params_out(:,5);
                Results.Sigma_SE=sqrt(CRLB_out(:,5));
            case 'XYZNB'
                Results.Z=Params_out(:,5);
                Results.Z_SE=sqrt(CRLB_out(:,5));
            case 'XYNBSXSY'
                Results.SigmaY=Params_out(:,5);
                Results.SigmaX=Params_out(:,6);
                Results.SigmaY_SE=sqrt(CRLB_out(:,5));
                Results.SigmaX_SE=sqrt(CRLB_out(:,6));
        end
        Results.LogLikelihood=LL_out;
        
        %Calculate p values
        Results.PValue=smi_core.GaussMLE.pValue(NP,size(obj.Data,1), ...
            Results.LogLikelihood);
        
        Statistics.FitTime=toc;
        Statistics.FitPerSecond=N/Statistics.FitTime;
        Statistics.NKernelsCalls = NKernelsCalls;
        Statistics.NFitsPerCall=NFitsPerCall;
        end
        
    end
    
    methods (Static)
        function [PValue]=pValue(NParams,FitBoxSize,LLR)
        %pValue Calculate p values from Log Likelihood Ratio
        %
        % INPUTS 
        %   NParams:    Number of fit parameters
        %   FitBoxSize: Linear size of fit box (Pixels)
        %   LLR:        Log Likelihood Ratio
        % OUTPUTS
        %   PValue:     P value of fit
        %
            X2_CDF=@(k,x)gammainc(x/2,k/2);
            K=FitBoxSize^2-NParams;
            X2=-2*LLR;
            PValue=1-X2_CDF(K,X2);
        end
    end
end
