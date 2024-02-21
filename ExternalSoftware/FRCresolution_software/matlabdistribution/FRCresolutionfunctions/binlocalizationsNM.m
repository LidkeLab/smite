%BINLOCALIZATIONSNM   Create image of binned localizations
%
% SYNOPSIS:
%   out = binlocalizationsNM(positions,xsize,ysize,pixelsize)
%
% PARAMETERS:
%   xsize
%      Size in x-direction of output (in nm)
%   ysize
%      Size in y-direction of output (in nm)
%   pixelsize
%      Pixel size of the output image (in nm)
%
% DEFAULTS:
%   xsize : image edge corresponds to maximum x-coordinate
%   ysize : image edge corresponds to maximum y-coordinate
%   pixelsize  = 10
%
% NOTES:
%  Columns of positions are assumed to be in nm.
%  Localizations with x or y-values smaller than -0.5*pixelsize are discarded.

% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen, Nov 2012


function out = binlocalizationsNM(varargin)

d = struct('menu','FRC resolution',...
           'display','Bin localizations in nm',...
           'inparams',struct('name',       {'positions',        'xsize',            'ysize',                        'ps'},...
                             'description',{'Localizations',    'X-size of output (nm)', 'Y-size of output (nm)',   'Superresolution pixel size (nm)'},...
                             'type',       {'array',            'array',            'array',                        'array'},...
                             'dim_check',  {{[],[-1 2],[-1 3]}, {[],0},             {[],0},                         {[],0}},...
                             'range_check',{'R',                'N+',            	'N+',                           [eps Inf]},...
                             'required',   {0,                  0,                  0,                              0},...
                             'default',    {[],                 [],                 [],                             10}...
                              ),...
           'outparams',struct('name',{'out'},...
                              'description',{'Image with binned localizations'},...
                              'type',{'array'}...
                              )...
           );       

if nargin == 1
   s = varargin{1};
   if ischar(s) & strcmp(s,'DIP_GetParamList')
      out = d;
      return
   end
end

try
   [positions,xsize,ysize,ps] = getparams(d,varargin{:});
catch
   if ~isempty(paramerror)
      error(paramerror)
   else
      error('Parameter parsing was unsuccessful.')
   end
end       

% Check that an output image size can be defined
if isempty(positions)
    if isempty(xsize) || isempty(ysize)
        error('No output image size can be defined.')
    else
        out = newim(xsize,ysize);
        return
    end
end

if isempty(ps)
    ps = 10;
end

zm = 1/ps;

if ~isempty(xsize)
    xsize = xsize/ps;
end

if ~isempty(xsize)
    ysize = ysize/ps;
end

out = binlocalizations(positions,xsize,ysize,zm);
