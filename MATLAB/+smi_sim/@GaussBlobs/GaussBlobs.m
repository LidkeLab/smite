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
%   smi_cuda_gaussBlobROIStack.ptx
%   smi_cuda_gaussBlobROIStack.cu


properties 
end

methods
end

methods (Static)
    [Model,Data] = gaussBlobROIStack(SZ,SMD,VarianceIm,Covariance,PixType)
    %[Model,Data] = gaussBlobImage(SZ,NFrames,SMD,Background,Density,VarianceIm)
    [Model,Data] = gaussBlobImage(SMD,SMF,Bg,Density)
    [BlobStack]=genRandomBlobImage(SZ,NFrames,Rho,Photons,PSFSigma,Bg)
    
    
    function success = unitTest()
    %unitTest Tests static methods using default parameters 
    % REQUIRES:
    %   DipImage
        fprintf('Testing gaussBlobROIStack...\n')
        success = 0;
        SaveDir = fullfile(tempdir, 'smite', 'unitTest', 'GaussBlobs');
        if ~isfolder(SaveDir)
           mkdir(fullfile(tempdir, 'smite'));
           mkdir(fullfile(tempdir, 'smite', 'unitTest'));
           mkdir(fullfile(tempdir, 'smite', 'unitTest', 'GaussBlobs'));
        end
        [Model,Data]=smi_sim.GaussBlobs.gaussBlobROIStack();
        figure; imagesc(sum(Model, 3)); colormap(gca, gray(256));
        saveas(gcf, fullfile(SaveDir, 'GB1.png'));
        figure; imagesc(sum(Data, 3)); colormap(gca, gray(256));
        saveas(gcf, fullfile(SaveDir, 'GB2.png'));
        success = 1;
        %dipshow(Model)
        %dipshow(Data)
%       fprintf('Testing gaussBlobImage...\n')
%       SMF = smi_core.SingleMoleculeFitting();
%       SMD = smi_core.SingleMoleculeData.createSMD();
%       [Model,Data]=smi_sim.GaussBlobs.gaussBlobImage(SMD, SMF);
%       figure; imagesc(sum(Model, 3)); colormap(gca, gray(256));
%       figure; imagesc(sum(Data, 3)); colormap(gca, gray(256));
        %dipshow(Model)
        %dipshow(Data)
    end
end
    
end
