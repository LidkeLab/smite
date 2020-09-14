function [resolution, frc_curve, resolution_hi, resolution_lo, ...
          resolution_t] = posToResolution(obj, SMR)
% posToResolution computes the image resolution from a list of localizations.
% The resolution is computed by taking blocks of localizations and assigning
% them randomly to half data sets.  The localizations in these sets are binned
% into 2 images from which the FRC curve and subsequently the resolution are
% computed.
%   
% INPUTS:
%    SMR             SMR structure with the following fields:
%       YSize        size (pixels) of the image
%       X,Y          (x, y) coordinates (pixels) of each localization [Nx1]
%       DatasetNum   dataset number of each localization              [Nx1]
%       FrameNum     frame number of each localization                [Nx1]
%       NFrames      number of frames in each dataset
%    obj:
%       SRImageZoom  image magnification factor
%       Repeats      number of times the FRC curve is computed for averaging
%
% OUTPUTS:
%    resolution      mean image resolution (pixels)
%    frc_curve       FRC curve vs. spatial frequency (1/nm)
%    resolution_hi   resolution + 1 standard deviation (pixels)
%    resolution_lo   resolution - 1 standard deviation (pixels)
%    resolution_t    mean of resolution over time
%
% REQUIRES:
%    DIPlib Image Resolution add-on
%    Parallel Processing Toolbox
%    NVidia GPU
%
% CITATION:
%    http://www.diplib.org/add-ons
%    Image Resolution, Reference: R.P.J. Nieuwenhuizen, K.A. Lidke, M. Bates,
%    D. Leyton Puig, D. Grünwald, S. Stallinga, B. Rieger, Measuring Image
%    Resolution in Optical Nanoscopy, Nature Methods, 10(6):557-562, 2013.
%    Also, see FRC_unitTest.m .

   % Compute absolute frame numbers as if there is just one dataset.
   FrameAbs = double(SMR.NFrames)*(SMR.DatasetNum - 1) + double(SMR.FrameNum);
   SRImageZoom = obj.SRImageZoom;;
   % Number of time blocks into which the dataset is split.  2 blocks is the
   % minimum accepted amount.
   blocks = 50;
   timefractions = 1;
   reps = obj.Repeats;

   [resolution, frc_curve, resolution_hi, resolution_lo, resolution_t] = ...
      postoresolution([SMR.X, SMR.Y, FrameAbs], SMR.YSize*SRImageZoom,   ...
                      SRImageZoom, timefractions, reps);

end
