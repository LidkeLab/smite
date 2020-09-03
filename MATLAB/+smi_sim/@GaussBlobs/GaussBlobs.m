 classdef GaussBlobs 
%GaussBlobs A collection of methods for generating 2D Gaussian Blob Images
% 
%   This class consists of two static methods that can be used to generate
%   image stacks containing Gaussian blobs. 
%   
%   gaussBlobROIStack generates a stack of images, each containing a 
%   single blob.
%   
%   gaussBlobImage generates a stack of images with multiple blobs and
%   internally uses gaussBlobROIStack.  This function is used in simulation
%   and display of 2D single molecule data. 
%
% REQUIRES:
%   MATLAB 2014a or later versions
%   Parallel Procesing Toolbox
%   NVidia GPU
%   cuda_gaussBlobROIStack.ptx
%   cuda_gaussBlobROIStack.cu


properties 
end

methods
end

methods (Static)
    [Model,Data] = gaussBlobROIStack(SZ,SMD,VarianceIm,Covariance,PixType)
    [Model,Data] = gaussBlobImage(SZ,NFrames,SMD,Background,Density,VarianceIm)
    
    function unitTest()
    %unitTest Tests static methods using default parameters 
    % REQUIRES:
    %   DipImage
        fprintf('Testing gaussBlobROIStack...\n')
        [Model,Data]=smi_sim.GaussBlobs.gaussBlobROIStack();
        dipshow(Model)
        dipshow(Data)
        fprintf('Testing gaussBlobImage...\n')
        [Model,Data]=smi_sim.GaussBlobs.gaussBlobImage();
        dipshow(Model)
        dipshow(Data)
    end
end
    
end