function frc_curve = posToFRC(obj, SMR)
% posToFRC calculates the Fourier Ring Correlation curve.
% The FRC curve is computed by taking blocks of localizations and assigning 
% them randomly to half data sets. The localizations in these sets are binned
% into 2 images from which the FRC curve is computed.
%   
% INPUTS:
%    SMR:            SMR structure with the following fields:
%       YSize        size (pixels) of the image
%       X,Y          (x, y) coordinates (pixels) of each localization [Nx1]
%       DatasetNum   dataset number of each localization              [Nx1]
%       FrameNum     frame number of each localization                [Nx1]
%       NFrames      number of frames in each dataset
%    obj:
%       SRImageZoom  image magnification factor
%
% OUTPUTS:
%    frc_curve       FRC curve vs. spatial frequency (1/nm)
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
   SRImageZoom = obj.SRImageZoom;
   % Number of time blocks into which the dataset is split.  2 blocks is the
   % minimum accepted amount.
   blocks = 50;

   frc_curve = postofrc([SMR.X, SMR.Y, FrameAbs], SMR.YSize*SRImageZoom, ...
                        SRImageZoom, blocks);

end
