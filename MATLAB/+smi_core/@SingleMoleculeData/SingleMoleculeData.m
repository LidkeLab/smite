classdef SingleMoleculeData
    % SingleMoleculeData A class defining the Single Molecule Data structure
    %
    % This datatype is one of the primary results structures in the smite
    % enviroment. The SMD structure is an input and output of many smi
    % methods. It intended to be extensible.
    % The SMD class implements tools for working with SMD structures,
    % but the data strcuture itself is not an object of the class.
    %
    % The structure has the following Properties:
    %
    % SMD:
    %   NDims:          Number of dimensions in localization information (2,3)
    %   NFrames:        Number of image frames in raw data sets
    %   NDatasets:      Number of 3D image stacks
    %   FrameRate:      Acquisition frame rate (1/seconds)
    %   PixelSize       Pixel size of camera projected onto sample (micrometers)
    %   XSize:          Number of pixels in X dimension of raw data
    %   YSize:          Number of pixels in Y dimension of raw data
    %   XBoxCorner:     X coordinate of top right box corner
    %   YBoxCorner:     Y coordinate of top right box corner
    %   ZOffset:        Z position of focal plane of sequence
    %   X:              Estimated X position
    %   Y:              Estimated Y position
    %   Z:              Estimated Z position
    %   Photons:        Estimated Photons  (Integrated collected photons)
    %   Bg:             Estimated Background (Photons/Pixel)
    %   PSFSigma:       Estimated or Fixed Sigma of 2D Gaussian PSF Model (symmetric PSF)
    %   PSFSigmaX:      Estimated or FixedX Sigma of 2D Gaussian PSF Model (asymmetric PSF)
    %   PSFSigmaY:      Estimated or FixedY Sigma of 2D Gaussian PSF Model (asymmetric PSF)
    %   X_SE:           Standard Error of X
    %   Y_SE:           Standard Error of Y
    %   Z_SE:           Standard Error of Z
    %   Photons_SE:     Standard Error of Photons
    %   Bg_SE:          Standard Error of Bg
    %   PSFSigma_SE:    Standard Error of PSFSigma
    %   PSFSigmaX_SE:   Standard Error of PSFSigmaX
    %   PSFSigmaY_SE:   Standard Error of PSFSigmaY
    %   DatasetNum:     File number from which localization originates
    %   FrameNum:       Frame number from which localization originates
    %   PValue:         p-value of fit
    %   LogLikelihood:  Log likelihood of fit
    %   ConnectID:      Identifies the same emitter accross multiple frames
    %   IndSMD:         Indices in original SMD corresponding to frame
    %                   connected localizations (e.g., indices in SMD
    %                   corresponding to localizations in SMDCombined).
    %   ThreshFlag:     Indicates a valid fit.  0=valid.  See SMA_Core.ThresholdSM
    %   DriftX:         X drift relative to first frame (Pixels) (NFrames x NDatasets)
    %   DriftY:         Y drift relative to first frame (Pixels) (NFrames x NDatasets)
    %   DriftZ:         Z drift relative to first frame (Pixels) (NFrames x NDatasets)
    %   IsTransformed:   Flag indicating channel reg. was performed on this SMD
    %
    % SEE ALSO:
    %   smi_core.SMF, smi_core.TR
    %
    
    properties
        
    end
    
    methods (Static)
        function [SMD] = createSMD()
            %createSMD Creates an empty Single Molecule Data (SMD) structure
            SMD.DatasetNum=[];
            SMD.FrameNum=[];
            SMD.X=[];
            SMD.Y=[];
            SMD.Z=[];
            SMD.X_SE=[];
            SMD.Y_SE=[];
            SMD.Z_SE=[];
            SMD.Photons=[];
            SMD.Photons_SE=[];
            SMD.Bg=[];
            SMD.Bg_SE=[];
            SMD.PSFSigma=[];
            SMD.PSFSigmaX=[];
            SMD.PSFSigmaY=[];
            SMD.PSFSigma_SE=[];
            SMD.PSFSigmaX_SE=[];
            SMD.PSFSigmaY_SE=[];
            SMD.PValue=[];
            SMD.LogLikelihood=[];
            SMD.ThreshFlag=[];
            SMD.ConnectID=uint32([]);
            SMD.IndSMD={};
            SMD.DriftX=[];
            SMD.DriftY=[];
            SMD.DriftZ=[];
            SMD.XBoxCorner=[];
            SMD.YBoxCorner=[];
            SMD.ZOffset=[];
            SMD.NDims=[];
            SMD.NFrames=[];
            SMD.NDatasets=[];
            SMD.FrameRate=[];
            SMD.PixelSize=[];
            SMD.XSize=[];
            SMD.YSize=[];
            SMD.IsTransformed=[];
        end
        
        [SMD] = catSMD(SMD1, SMD2)
        [Density] = computeDensity(SMD);
        [SMD] = isolateSubSMD(SMD, SubIndices);
        
    end
end