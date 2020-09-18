%FPC_SLICES   Compute the FPC from 2 images in 3 planes
%
% First a Tukey window is applied to the 2 images, and subsequently the FPC
% is computed from the Fourier Transforms of the two images, for slices
% where the x,y, or z-component of the frequency vector is 0.
%
% SYNOPSIS:
%   [fpc_xy fpc_xz fpc_yz] = fpc_slices(in1,in2)
%
% NOTES:
%   Non-cubed images are zero padded to make them square.
%   Images with even sizes are zero padded to make the image sizes uneven.

% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen, Oct 2012

function varargout = fpc_slices(varargin)

d = struct('menu','FRC resolution',...
           'display','FLC',...
           'inparams',struct('name',       {'in1',    'in2'},...
                             'description',{'Image 1','Image 2'},...
                             'type',       {'image',  'image'},...
                             'dim_check',  {3,        3},...
                             'range_check',{[],       []},...
                             'required',   {1,        1},...
                             'default',    {'in1',    'in2'}...
                              ),...
           'outparams',struct('name',{'fpc_xy','fpc_xz','fpc_yz'},...
                              'description',{'FPC xy-slice','FPC xz-slice','FPC yz-slice'},...
                              'type',{'image','image','image'}...
                              )...
           );       
       
if nargin == 1
   s = varargin{1};
   if ischar(s) && strcmp(s,'DIP_GetParamList')
      varargout{1} = struct('menu','none');
      return
   end
end

try
   [in1, in2] = getparams(d,varargin{:});
catch
   if ~isempty(paramerror)
      error(paramerror)
   else
      error('Parameter parsing was unsuccessful.')
   end
end

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

% Fourier transform input images
f1 = ft(in1);
f2 = ft(in2);
clear in1 in2

f1 = f1*(rr(sz)<=floor(sz(1)/2));
f2 = f2*(rr(sz)<=floor(sz(1)/2));

% Compute correlations
h12 = real(f1*conj(f2)); % the product is complex but later the complex part cancels
                         % out in the sum over the planes, faster to due it this way.
h11 = real(f1*conj(f1)); % just change the datatype as output is real
h22 = real(f2*conj(f2));
clear f1 f2

% Calculate FPC in the qx-qy plane by summing the z dimension and taking line sums
fpc_num = wkk(squeeze(sum(h12,[],3)));                                          % Numerator of the FPC in the qx-qy plane
fpc_denom = sqrt(wkk(squeeze(sum(h11,[],3)))*wkk(squeeze(sum(h22,[],3))));      % Denominator of the FPC in the qx-qy plane
fpc_xy = real(fpc_num./fpc_denom);

% Calculate FPC in the qx-qz plane
fpc_num = wkk(squeeze(sum(h12,[],2)));                                          % Numerator of the FPC in the qx-qz plane
fpc_denom = sqrt(wkk(squeeze(sum(h11,[],2)))*wkk(squeeze(sum(h22,[],2))));      % Denominator of the FPC in the qx-qz plane
fpc_xz = real(fpc_num./fpc_denom);

% Calculate FPC in the qy-qz plane
fpc_num = wkk(squeeze(sum(h12,[],1)));                                          % Numerator of the FPC in the qy-qz plane
fpc_denom = sqrt(wkk(squeeze(sum(h11,[],1)))*wkk(squeeze(sum(h22,[],1))));      % Denominator of the FPC in the qy-qz plane
fpc_yz = real(fpc_num./fpc_denom);

varargout{1} = fpc_xy;
varargout{2} = fpc_xz;
varargout{3} = fpc_yz;

