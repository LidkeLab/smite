%FRCRESOLUTION   Compute the resolution from 2 images
%
% SYNOPSIS:
%   [resolution frc_out resolution_high resolution_low] = frcresolution(in1,in2)
%
% NOTES:
%   Non-square images are zero padded to make them square.

% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen, Oct 2012

function varargout = frcresolution(varargin)

d = struct('menu','FRC resolution',...
           'display','Resolution from images',...
           'inparams',struct('name',       {'in1',      'in2',},...
                             'description',{'Image 1',  'Image 2',},...
                             'type',       {'image',    'image',},...
                             'dim_check',  {2,          2,},...
                             'range_check',{[],         [],},...
                             'required',   {1,          1,},...
                             'default',    {[],         [],}...
                              ),...
           'outparams',struct('name',{'resolution','frc_out','resolution_high','resolution_low'},...
                              'description',{'Resolution','FRC curve','Upper bound','Lower bound'},...
                              'type',{'array','array','array','array'}...
                              )...
           );       
       
if nargin == 1
   s = varargin{1};
   if ischar(s) & strcmp(s,'DIP_GetParamList')
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

% Calculate FRC curve
frc_out = frc(in1,in2);

% Calculate the resolution
sz = max(imsize(in1));

if (~sum(in1)) || (~sum(in2))
    varargout{1} = NaN;
    varargout{2} = radialsum(0*in1);
    varargout{3} = NaN;
    varargout{4} = NaN;
    return;
end

switch nargout
    case 1
        resolution = frctoresolution(frc_out,sz);
        varargout{1} = resolution;
        return
    case 2
        resolution = frctoresolution(frc_out,sz);
        varargout{1} = resolution;
        varargout{2} = frc_out;
        return        
    case 3
        [resolution resolution_high] = frctoresolution(frc_out,sz);
        varargout{1} = resolution;
        varargout{2} = frc_out;
        varargout{3} = resolution_high;
        return
    case 4
        [resolution resolution_high resolution_low] = frctoresolution(frc_out,sz);
        varargout{1} = resolution;
        varargout{2} = frc_out;
        varargout{3} = resolution_high;
        varargout{4} = resolution_low;
        return
end
