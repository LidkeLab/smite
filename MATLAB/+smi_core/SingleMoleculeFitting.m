classdef SingleMoleculeFitting<handle
% SingleMoleculeFitting A class defining the Single Molecule Fitting structure
%
% The SMF structure is a structure of structures that collectively contain
% all parameters required to go from raw data to an SMD results structure.
% The SMF structure is an input of many smi methods. It
% intended to be extensible to enable new analysis tools and methods.
% The SMF class implements tools for working with SMF structures,
% but the data structure itself is not an object of the class.
%
% Parameters of sub-structures are explained in more detail in
% the classes and methods that use them.  An incomplete list of classes
% that use each sub-structure is listed in {}.
%
% The SMF structure has the following sub-structures and fields:
%
% SMF:  Fields that are structures: 
%
% Data:             {LoadData}
%   FileName:       File name (char array or cell array of char array)
%   FileDir:        File directory (char array or cell array of char array)
%   ResultsDir:     Results directory (char array)(Default='FileDir/Results')
%   CameraType:     'EMCCD','SCMOS' (Default='EMCCD')
%   CameraGain:     Camera Gain, scalar or image (Default=1)
%   CameraOffset:   Camera Offset, scalar or image (Default=0)
%   CameraNoise:    Camera readnoise, scalar or image (Default=0)
%   FrameRate:      Data Collection Frame Rate (1/s) 
%   PixelSize:      Camera back-projected pixel size (micron)   
%
% BoxFinding:       {FindROI}
%   BoxSize:        Linear box size for fitting (Pixels)(Default=7)
%   BoxOverlap:     Overlap of boxes allowed (Pixels)(Default=2)
%   MinPhotons:     Minimum number of photons from emitter (Default=200)
%
% Fitting           {GaussMLE}
%   PSFSigma:   Initial or fixed Sigma of 2D Gaussian PSF Model (Pixels)(Default=1)
%   FitType:    See fit class for options  (Default='XYNB')
%   Iterations: Newton Raphson iterations (Default=20)
%   ZFitStruct: Structure for astigmatic fitting:
%       Ax:         Asigmatism fit parameter (see GaussMLE)         
%       Ay:         Asigmatism fit parameter (see GaussMLE)
%       Bx:         Asigmatism fit parameter (see GaussMLE)
%       By:         Asigmatism fit parameter (see GaussMLE)
%       Gamma:      Asigmatism fit parameter (see GaussMLE)
%       D:          Asigmatism fit parameter (see GaussMLE)
%
% Thresholding      {ThresholdFits,SRA}
%   MaxSE_XY:       Maximum allowed precision in x,y (Pixels)(Default=.2)
%   MaxZ_SE:        Maximum allowed precision in z (Microns)(Default=.5)
%   LoS:            Minimum accepted p-value from fit (Default=.01)
%   MinPSFSigma:    Minimum PSF Sigma from fit (Pixels)(Default=.5);
%   MaxPSFSigma:    Maximum PSF Sigma from fit (Pixels)(Default=2);
%   MinPhotons:     Minimum accepted photons from fit (Default=100)
%   MaxBg:          Maximum background accepted from fit (Default=Inf)
%
% Threshold         {Threshold,SRA}
%   MinMax.X:       Min/Max x-position (Pixels)(Default=[])
%   MinMax.Y:       Min/Max y-position (Pixels)(Default=[])
%   MinMax.Z:       Min/Max z-position (Pixels)(Default=[])
%   MinMax.Photons: Min/Max photons accepted from fit (Default=[100,Inf])
%   MinMax.Bg:      Min/Max background accepted from fit (Default=[])
%   MinMax.PSFSigma:Min/Max PSF Sigma from fit (Pixels)(Default=[.5,2])
%   MinMax.X_SE:    Min/Max allowed precision in x (Pixels)(Default[0,0.2])
%   MinMax.Y_SE:    Min/Max allowed precision in y (Pixels)(Default=[0,0.2])
%   MinMax.Z_SE:    Min/Max allowed precision in z (Microns)(Default=[0,0.5])
%   MinMax.Pvalue:  Min/Max accepted p-value from fit (Default=[.01,1])
%   MeanPhotons:    Mean photons (Default=200) (see FindROI)
%
% FrameConnection:  {FrameConnect,SRA}
%   MaxSeparation:  Maximum separation for connection (Pixels)(Default=1)
%   MaxFrameGap:    Maximum frame gap for connection (Frames)(Default=4)
%   LoS:            Minimum accepted p-value for connection (Default=.01)
%
% DriftCorrection   {DriftCorrection}
% 
% Tracking          {SPT}
%   Method:         Type of method used for tracking (Default='smi_spt')
%   D:              Diffusion Constant (Pixels^2/Frame) (Default=1)
%   K_on:           Off to On Rate (Frame^-1) (Default=.1)
%   K_off:          On to Off Rate (Frame^-1) (Default=.1)
%   MaxFrameGap:    Maximum frame gap for Gap Closing (Pixels) (Default=10)
%   MaxDist:        Maximum distance gap for Gap Closing (Pixels) (Default=10)
%   MinTrackLength  Minimum track length of trajectory (Frames) (Default=3)
%

% created by:
% Keith Lidke, Hanieh Mazloom-Farsibaf, David Schodt. Lidke Lab 2018

    
    
    methods (Static)
        
        function [SMF] = createSMF()
            %createSMF Creates an default Single Molecule Fitting (SMF) structure
            %
            % INPUTS:
            %   (none)
            % OUTPUTS:
            %   SMF:    An SMF structure with all fields set to defaults
            % REQUIRES:
            %   (none)
            %
            
            SMF.Data.CameraType='EMCCD';
            SMF.Data.CameraGain=1;
            SMF.Data.CameraOffset=0;
            SMF.Data.CameraReadNoise=0;
            
            SMF.BoxFinding.BoxSize=7;
            SMF.BoxFinding.BoxOverlap=2;
            SMF.BoxFinding.MinPhotons=200;
            
            %Fitting
            SMF.Fitting.PSFSigma=1;
            SMF.Fitting.FitType='XYNB';
            SMF.Fitting.Iterations=20;
            SMF.Fitting.ZFitStruct.Ax=[];
            SMF.Fitting.ZFitStruct.Ay=[];
            SMF.Fitting.ZFitStruct.Bx=[];
            SMF.Fitting.ZFitStruct.By=[];
            SMF.Fitting.ZFitStruct.Gamma=[];
            SMF.Fitting.ZFitStruct.D=[];

            %Thresholding
            SMF.Thresholding.MaxXY_SE=.2;
            SMF.Thresholding.MaxZ_SE=.05;
            SMF.Thresholding.MinPValue=.01;
            SMF.Thresholding.MinPSFSigma=0.5;
            SMF.Thresholding.MaxPSFSigma=2;
            SMF.Thresholding.MeanPhotons=200;
            SMF.Thresholding.MinPhotons=100;
            SMF.Thresholding.MaxBG=Inf;

            %Threshold
            SMF.Threshold.MinMax.X=[];
            SMF.Threshold.MinMax.Y=[];
            SMF.Threshold.MinMax.Z=[];
            SMF.Threshold.MinMax.Photons=[100, Inf];
            SMF.Threshold.MinMax.Bg=[];
            SMF.Threshold.MinMax.PSFSigma=[0.5, 2];
            SMF.Threshold.MinMax.X_SE=[0, 0.2];
            SMF.Threshold.MinMax.Y_SE=[0, 0.2];
            SMF.Threshold.MinMax.Z_SE=[0, 0.05];
            SMF.Threshold.MinMax.Pvalue=[0.01, 1];
            SMF.Threshold.MeanPhotons=200;
            
            %FrameConnection
            SMF.FrameConnection.MaxSeparation=1; % pixels 
            SMF.FrameConnection.MaxFrameGap=4; % frames
            SMF.FrameConnection.LoS=.01;

            %DriftCorrection
            
            %Tracking
            SMF.Tracking.TrackMethods='SMA_SPT';
            SMF.Tracking.D=1;
            SMF.Tracking.K_on=.1;
            SMF.Tracking.K_off=.1;
            SMF.Tracking.MaxFrameGap=10;
            SMF.Tracking.MaxDist=10;
            SMF.Tracking.MinTrackLength=3;
        end
    end
end
