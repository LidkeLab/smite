%FPC   Compute FPC from 2 images
%
% First a Tukey window is applied to the 2 images, and subsequently the FPC
% is computed from the Fourier Transforms of the two images.
%
% SYNOPSIS:
%   [fpc_xy fpc_xz fpc_yz fpc_out] = fpc(in1,in2,angles,maxfreq)
%
%   angles
%      Number of different orienations of the orientation planes. 
%   maxfreq
%      A cubic volume in Fourier space from -maxfreq to maxfreq is cut out
%      for the FPC calculation
%
% DEFAULTS:
%   angles = 64
%   maxfreq = 0.5
%
% NOTES:
%   Non-cubed images are zero padded to make them square.
%   Images with even sizes are zero padded to make the image sizes uneven.
%   The variable angles is rounded down to the nearest power of 4.

% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Bernd Rieger and Robert Nieuwenhuizen, Oct 2012

function varargout = fpc(varargin)

d = struct('menu','FRC resolution',...
           'display','FPC',...
           'inparams',struct('name',       {'in1',    'in2',        'N_angles',         'maxfreq',              'full_fpc'},...
                             'description',{'Image 1','Image 2',    'Rotation angles',  'Maximum frequency',    'Compute full FPC'},...
                             'type',       {'image',  'image',      'array',            'array',                'boolean'},...
                             'dim_check',  {3,        3,            {0},                {0},                    0},...
                             'range_check',{[],       [],           'N+',               [eps 0.5],              []},...
                             'required',   {1,        1,            0,                  0,                      0},...
                             'default',    {'in1',    'in2',        64,                 0.5,                    0}...
                              ),...
           'outparams',struct('name',{'fpc_xy','fpc_xz','fpc_yz','fpc_out'},...
                              'description',{'FPC in xy-plane','FPC in xz-plane','FPC in yz-plane','Full FPC'},...
                              'type',{'image','image','image','image'}...
                              )...
           );       
       
if nargin == 1
   s = varargin{1};
   if ischar(s) && strcmp(s,'DIP_GetParamList')
      varargout{1} = d;
      return
   end
end

try
   [in1, in2, N_angles, maxfreq, full_fpc] = getparams(d,varargin{:});
catch
   if ~isempty(paramerror)
      error(paramerror)
   else
      error('Parameter parsing was unsuccessful.')
   end
end

% Only compute slices of the FPC if 3 output arguments are used
if nargout>4
    error('Incorrect number of output arguments.')
end

if nargout>1
    [fpc_xy fpc_xz fpc_yz] = fpc_slices(in1,in2);
    varargout{1} = fpc_xy;
    varargout{2} = fpc_xz;
    varargout{3} = fpc_yz;
    if nargout<4
        return
    end
    
    if ~full_fpc
        varargout{4} = [];
        return
    end
end

% Check that N_angles exists and is 4 or larger
if N_angles<4
    N_angles = 4;
end

% Check that maxfreq exists
sz = imsize(in1);

if any(imsize(in2) ~= sz)
    error('Image 1 and image 2 have different image sizes.');
end

%% Mask and extend images

% Compute mask in x-direction
nfac = 8;                                                   % Image width / Width of edge region
x_im = xx(sz(1),sz(2),sz(3))/sz(1);
mask = 0.5-0.5*cos(pi*nfac*x_im);          
mask(abs(x_im)<((nfac-2)/(nfac*2))) = 1;
clear x_im

% Check that input images are square and mask
if sz(1) == sz(2)
    mask = mask*rot90(mask);
else
    % Compute mask in y-direction
    y_im = yy(sz(1),sz(2))/sz(2);
    mask_y = 0.5-0.5*cos(pi*nfac*y_im);          
    mask_y(abs(y_im)<((nfac-2)/(nfac*2))) = 1;
    mask = mask*mask_y;
    clear mask_y
end

% Mask input images
in1 = mask*in1;
in2 = mask*in2;
clear mask

% Make images cubed through zero padding if necessary
sz(sz<max(sz)) = max(sz);
in1 = extend(in1,sz);
in2 = extend(in2,sz);
    
% Make image size uneven through zero padding
if ~mod(sz(1),2)
    sz = sz+1;
    in1 = extend(in1,sz);
    in2 = extend(in2,sz);
end

%% Evaluate Fourier Plane Correlation

% Get all plane orientations (polar angle and azimuthal angle)
[pol_ang,azim_ang] = sphere_tesselation(floor(log(N_angles)/log(4)));

% Find the width in Fourier space of the pixels that are used
n = ceil(2*maxfreq*sz(1));                              % 2 times the radius in Fourier space of highest unmasked frequency component
if (~mod(n,2) && n>1)
    n = n-1;                                            % Make n uneven if necessary
end

% Fourier Plane Correlation
fpc_raw = gfca3D_sub(in1,in2,pol_ang,azim_ang,n);       % FPC for the specified angles and radii in Fourier space
fpc_raw = im2mat(fpc_raw)';

% Find xyz-coordinates of plane correlation values
Np = size(fpc_raw,1);                                   % Number of Fourier plane radii
qradius = ceil((0:Np-1)-Np/2)';                         % Plane radii
[~, pol_ang] = ndgrid(qradius,pol_ang);            
[q_mat, azim_ang] = ndgrid(qradius,azim_ang);
[x,y,z] = sph2cart(azim_ang,pi/2-pol_ang,q_mat);
clear q_mat pol_ang azim_ang

% Define grid of output values
[xc,yc,zc] = meshgrid(qradius,qradius,qradius); 

% Check if a warning will be displayed when calling TriScatteredInterp
s = warning('query', 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId');
if strcmp(s.state,'on')
    warning('off', 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId');
end

% Interpolate the FPC data on the gridpoints (xc,yc,zc)
x = reshape(x,[numel(x) 1]);
y = reshape(y,[numel(y) 1]);
z = reshape(z,[numel(z) 1]);
fpc_raw = reshape(fpc_raw,[numel(fpc_raw) 1]);
Int_obj = TriScatteredInterp(x,y,z,fpc_raw);            % Interpolated hypersurface v = fpc_raw(x,y,z)
fpc_out = mat2im(Int_obj(xc,yc,zc));                    % Evaluate value of hypersurface at (xc,yc,zc)
fpc_out(isnan(fpc_out)) = 0;

% Return the state of the warning from griddata to its original state
if strcmp(s.state,'on')
    warning('on', 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId');
end

if nargout == 1
    varargout{1} = fpc_out;
else
    varargout{4} = fpc_out;
end