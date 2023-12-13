%FLC   Compute FLC from 2 images
%
% First a Tukey window is applied to the 2 images, and subsequently the FLC
% is computed from the Fourier Transforms of the two images. The
% maximum frequency for which the FLC is computed is 1/2.
%
% SYNOPSIS:
%   flc_out = flc(in1,in2)
%
% NOTES:
%   Non-square images are zero padded to make them square.
%   Images with even sizes are zero padded to make the image sizes uneven.

% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen, Oct 2012

function flc_out = flc(varargin)

d = struct('menu','FRC resolution',...
           'display','FLC',...
           'inparams',struct('name',       {'in1',    'in2'},...
                             'description',{'Image 1','Image 2'},...
                             'type',       {'image',  'image'},...
                             'dim_check',  {2,        2},...
                             'range_check',{[],       []},...
                             'required',   {1,        1},...
                             'default',    {'in1',    'in2'}...
                              ),...
           'outparams',struct('name',{'flc_out'},...
                              'description',{'FLC'},...
                              'type',{'image'}...
                              )...
           );       
       
if nargin == 1
   s = varargin{1};
   if ischar(s) && strcmp(s,'DIP_GetParamList')
      flc_out = d;
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

% Compute mask in x-direction
nfac = 8;                                                   % Image width / Width of edge region
x_im = xx(sz(1),sz(2))/sz(1);
mask = 0.5-0.5*cos(pi*nfac*x_im);          
mask(abs(x_im)<((nfac-2)/(nfac*2))) = 1;
clear x_im

% Check that input images are square and mask
if sz(1) == sz(2)
    mask = mask*rot90(mask);
    
    % Mask input images
    in1 = mask*in1;
    in2 = mask*in2;
else
    % Compute mask in y-direction
    y_im = yy(sz(1),sz(2))/sz(2);
    mask_y = 0.5-0.5*cos(pi*nfac*y_im);          
    mask_y(abs(y_im)<((nfac-2)/(nfac*2))) = 1;
    mask = mask*mask_y;
    clear mask_y
    
    % Mask input images
    in1 = mask*in1;
    in2 = mask*in2;
    
    % Make images square through zero padding
    sz(sz<max(sz)) = max(sz);
    in1 = extend(in1,sz);
    in2 = extend(in2,sz);
end
clear mask

% Make image size uneven through zero padding
if ~mod(sz(1),2)
    sz = sz+1;
    in1 = extend(in1,sz);
    in2 = extend(in2,sz);
end

% Fourier transform input images
in1 = ft(in1);
in2 = ft(in2);

% Mask frequencies larger than 0.5
in1 = in1*(rr(sz(1),sz(2))<=floor((sz(1)-1)/2));
in2 = in2*(rr(sz(1),sz(2))<=floor((sz(1)-1)/2));

% Compute fourier line correlation
flc_num = wkk(in1.*conj(in2));          % Numerator of the FLC
in1 = abs(in1).^2;
in2 = abs(in2).^2;
flc_denom = sqrt(wkk(in1)*wkk(in2));    % Denominator of the FLC
flc_out = real(flc_num./flc_denom);
