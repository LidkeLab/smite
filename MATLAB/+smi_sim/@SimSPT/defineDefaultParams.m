function [ParamStruct] = defineDefaultParams()
%defineDefaultParams creates a ParamStruct with all default values set.
% This method creates a structure of parameters 'ParamStruct' with all
% values set to some meaningful default value.
%
% OUTPUTS:
%   ParamStruct: structure array of default parameters.
%                ParticleDensity: Density of particles within the 
%                                 simulation. (particles / pixel^2)
%                                 (Default = 0.005)
%                NFrames: Number of (full) frames in the simulation (i.e.,
%                         if NFrames = 100 and SubframeDensity = 2, we will
%                         simulate 200 subframes). (Default = 100)
%                SubframeDensity: Number of subframes per frame (i.e., each 
%                                 frame of the simulation represents the
%                                 motion blur of the subframes). If 
%                                 mod(NFrames, NSubframes)~=0, remainder 
%                                 frames are thrown away. 
%                                 (1/frame)(Default = 1)
%                FrameSize: Size of the frame within which the trajectories 
%                           are simulated. (pixels)(Default = [32, 32])
%                PSFSigma: Standard deviation of the best Gaussian fit to
%                          simulated emitters. (pixels)(Default = 1.3)
%                InitialDensityMask: Binary mask defining the allowable 
%                                    region for initial particle placement.
%                                    (Default = ones(FrameSize))
%                LabelingEfficiency: Probability of a target being labeled.
%                                    (Default = 1)
%                BoundaryCondition: Boundary condition applied when 
%                                   particles reach a boundary specified by
%                                   FrameSize. This can be 'Periodic', 
%                                   'Reflecting', or 'Free'. 
%                                   (Default = 'Periodic')
%                Intensity: Photons present in the simulated trajectory in
%                           each full frame. 
%                           (photons / frame)(Default = 1000)
%                Bg: Uniform background signal per full frame stored for
%                    trajectories. (photons / frame)(Default = 5)
%                D: Diffusion coefficients(s) for the Brownian motion 
%                   trajectories. If the size of this array is greater than
%                   or equal to the number of trajectories, the 1:NTraj 
%                   trajectories will have a diffusion constant taken in 
%                   order from this array.  If this array is smaller than 
%                   the number of trajectories, the diffusion constant for
%                   each trajectory will be randomly sampled from this
%                   array. (pixels^2 / frame)(Default = 0.1)
%                InteractionDistance: Distance between particles in a
%                                     dimer.  Distances for oligomers may 
%                                     not be equal to this value.
%                                     (pixels)(Default = 0.5)
%                InteractionProb: Probability that two particles within 
%                                 InteractionDistance of each other in a 
%                                 given frame will be dimerized.  Note that
%                                 this probability is "tested" at each 
%                                 sub-frame of the simulation.
%                                 (Default = 0.5)
%                RestrictToDimers: Restrict oligomerization to only order 2
%                                  (dimers). (Default = true)
%                KDisconnect: Rate parameter for the disconnection of 
%                             oligomers. (1 / frame)(Default = 0.1)
%                KOnToBleach: Rate parameter defining photobleaching.
%                            (1 / frame)(Default = 0.1)
%                KOnToOff: Rate parameter defining the turning off of 
%                          visible emitters. (1 / frame)(Default = 0.2)
%                KOffToOn: Rate parameter defining the turning on of 
%                          dark emitters. (1 / frame)(Default = 0.8)
%                PMiss: Probability of missing a localization of a visible
%                       emitter. (Default = 0.01)
%                Bg: Bg specifies a uniform background present in the raw 
%                    data. (photons)(Default = 5)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Populate the output structure.
ParamStruct.ParticleDensity = 0.005;
ParamStruct.NFrames = 100;
ParamStruct.SubframeDensity = 1;
ParamStruct.FrameSize = [32, 32];
ParamStruct.PSFSigma = 1.3;
ParamStruct.InitialDensityMask = ones(ParamStruct.FrameSize);
ParamStruct.LabelingEfficiency = 1;
ParamStruct.BoundaryCondition = 'Periodic';
ParamStruct.Intensity = 1000;
ParamStruct.Bg = 5;
ParamStruct.D = 0.1;
ParamStruct.InteractionDistance = 0.5;
ParamStruct.InteractionProb = 0.5;
ParamStruct.RestrictToDimers = true;
ParamStruct.KDisconnect = 0.1;
ParamStruct.KOnToBleach = 0.01;
ParamStruct.KOnToOff = 0.2;
ParamStruct.KOffToOn = 0.8;
ParamStruct.PMiss = 0.01;


end