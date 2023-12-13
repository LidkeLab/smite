%GFCA3D_SUB   Compute FPC from 2 images for specified angles
%
% Computes the FPC for the given combination of angles and radii in Fourier
% space.
%
% SYNOPSIS:
%   fpc_out = gfca3D_sub(in1,in2,polar_angles,azim_angles,N_freq)
%
%   polar_angles
%      List of M polar angles of the planes' orientations. 
%   azim_angles
%      List of M azimuthal angles of the planes' orientations. 
%   N_freq
%      Width in pixels of the volume in frequency space that is used.
%
% NOTES:
%   Input images must be cubic.

% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Bernd Rieger and Robert Nieuwenhuizen, Oct 2012

function fpc_out = gfca3D_sub(varargin)

d = struct('menu','FRC resolution',...
           'display','FLC',...
           'inparams',struct('name',       {'in1',    'in2',        'pol_ang',           'azim_ang',        'Nsub'},...
                             'description',{'Image 1','Image 2',    'Polar angles',     'Azimuthal angles', 'Cutout volume'},...
                             'type',       {'image',  'image',      'array',            'array',            'array'},...
                             'dim_check',  {3,        3,            [-1 1],             [-1 1],             {[],0}},...
                             'range_check',{[],       [],           [],                 [],                 'N+'},...
                             'required',   {1,        1,            1,                  1,                  0},...
                             'default',    {'in1',    'in2',        [],                 [],                 []}...
                              ),...
           'outparams',struct('name',{'fpc_out'},...
                              'description',{'FPC'},...
                              'type',{'image'}...
                              )...
           );       
       
if nargin == 1
   s = varargin{1};
   if ischar(s) && strcmp(s,'DIP_GetParamList')
      fpc_out = struct('menu','none');
      return
   end
end

try
   [in1, in2, pol_ang, azim_ang, Nsub] = getparams(d,varargin{:});
catch
   if ~isempty(paramerror)
      error(paramerror)
   else
      error('Parameter parsing was unsuccessful.')
   end
end

%... it is assumed that in1 and in2 are NpxNpxNp images...
%... with isotropic sampling...
%... with Np preferably an odd number to guarantee symmetry in end result... 
sz = imsize(in1);

if any(imsize(in2) ~= sz)
    error('Image 1 and image 2 have different image sizes.');
end

if any(sz ~= max(sz))
    error('Input images are not cubic.');
end

% Check that the numbers of angles are the same
if any(size(pol_ang) ~= size(azim_ang))
    error('The number of polar and azimuthal angles must be the same.');
end

% Define the volume in frequency space that is used if necessary
if isempty(Nsub)
    Nsub = sz(1);
end

% Fourier transform input images
fprintf('Computing Fourier transformation\n')
f1 = ft(in1);
f2 = ft(in2);
clear in1 in2

% Cutout volume in Fourier space
f1 = cut(f1,Nsub);
f2 = cut(f2,Nsub);
Np = imsize(f1,1);
Ncenter = floor(Np/2);

% Compute correlations
h12 = real(f1*conj(f2)); % the product is complex but later the complex part cancels
                         % out in the sum over the planes, faster to due it this way.
h11 = real(f1*conj(f1)); % just change the datatype as output is real
h22 = real(f2*conj(f2));
clear f1 f2

% Prepare computation of rotation matrices for all directions to +e_z
cospol = cos(pol_ang);
sinpol = sin(pol_ang);
cosazi = cos(azim_ang);
sinazi = sin(azim_ang);

% Loop over all directions, for each direction rotation to +e_z and sum
% over rotated xy-plane
fprintf('Rotating..\n')
fpc_out = zeros(numel(pol_ang),Np);

matlabpool(4)
parfor jj=1:numel(pol_ang)
  fprintf('%d of %d\n',jj,numel(pol_ang))
  Rz = [cosazi(jj),sinazi(jj),0;-sinazi(jj),cosazi(jj),0;0,0,1];
  Ry = [cospol(jj),0,-sinpol(jj);0,1,0;sinpol(jj),0,cospol(jj)];
  R = Ry*Rz;
  
  % ...somehow we need the transpose of R here...
  h12rot = rotation3d(h12,'direct',R',[Ncenter, Ncenter, Ncenter]);
  h11rot = rotation3d(h11,'direct',R',[Ncenter, Ncenter, Ncenter]);
  h22rot = rotation3d(h22,'direct',R',[Ncenter, Ncenter, Ncenter]);
  
  h12av = squeeze(sum(h12rot,[],[1,2]));
  h11av = squeeze(sum(h11rot,[],[1,2]));
  h22av = squeeze(sum(h22rot,[],[1,2]));
  fpc_out(jj,:) = im2mat(h12av./sqrt(h11av)./sqrt(h22av));
end
matlabpool close
% Correct for q = 0
% fpc(:,Np/2+1) = 1.0;

fpc_out = mat2im(fpc_out);


