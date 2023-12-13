%BINLOCALIZATIONS   Create image of binned localizations
%
% SYNOPSIS:
%   out = binlocalizations(positions,xsize,ysize,zoomfactor)
%
% PARAMETERS:
%   xsize
%      Size in x-direction of output
%   ysize
%      Size in y-direction of output
%   zoomfactor
%      Pixels per unit distance in the positions list
%
% DEFAULTS:
%   xsize : image edge corresponds to maximum x-coordinate
%   ysize : image edge corresponds to maximum y-coordinate
%   zoomfactor  = 1
%
% NOTES:
%  Columns of positions are assumed to have the same dimensionality.
%  Localizations with x or y-values smaller than -0.5 are discarded.

% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen, Oct 2012


function out = binlocalizations(varargin)

d = struct('menu','FRC resolution',...
           'display','Bin localizations',...
           'inparams',struct('name',       {'positions',        'xsize',            'ysize',            'zm'},...
                             'description',{'Localizations',    'X-size of output', 'Y-size of output', 'Zoom'},...
                             'type',       {'array',            'array',            'array',            'array'},...
                             'dim_check',  {{[],[-1 2],[-1 3]}, {[],0},             {[],0},             {[],0}},...
                             'range_check',{'R',                'N+',            	'N+',               [eps Inf]},...
                             'required',   {0,                  0,                  0,                  0},...
                             'default',    {[],                 [],                 [],                 1}...
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
   [positions,xsize,ysize,zm] = getparams(d,varargin{:});
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

if isempty(zm)
    zm = 1;
end

% Offset and rescale positions
positions = single(positions);
positions(:,1:2) = positions(:,1:2) + 0.5;
positions(:,1:2) = positions(:,1:2)*zm;

% Define output sizes if necessary
if isempty(xsize)
    xsize = ceil(max(positions(:,1)));
end

if isempty(ysize)
    ysize = ceil(max(positions(:,2)));
end

% Filter localizations outside the FOV
keep = positions(:,1)>=0;
keep = keep & positions(:,2)>=0;
keep = keep & positions(:,1)<=xsize;
keep = keep & positions(:,2)<=ysize;
positions = positions(keep,:); 

% Bin localizations
out = cHistRecon(xsize,ysize,positions(:,1),positions(:,2),0)';
out = mat2im(out);
