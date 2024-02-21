%POSTOFRC   Compute FRC curve from a list of localizations
%
% The FRC curve is computed by taking blocks of localizations and assigning 
% them randomly to half data sets. The localizations in these sets are 
% binned into 2 images from which the FRC curve is computed.
%
% SYNOPSIS:
%   frc_out = postofrc(positions,size,zoomfactor,blocks)
%
%   size
%      Size of the binned images
%   zoom
%      Pixels per unit distance in the positions list in the binned images
%   blocks
%      Number of time blocks into which the dataset is split. 2 block is
%      the minimum accepted amount.
%
% DEFAULTS:
%   size : image edge corresponds to maximum y-coordinate
%   zoomfactor  = 1
%   blocks = 50
%
% NOTES:
%   The third column of positions is assumed to contain time stamps of the
%   localizations. If positions only has 2 columns then the localizations 
%   are assumed to be time-ordered.
%   The minimum number of time blocks is 2.

% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen, Oct 2012

function frc_out = postofrc(varargin)

d = struct('menu','FRC resolution',...
           'display','FRC from localizations',...
           'inparams',struct('name',       {'positions',        'sz',           'zm',           'nblocks'},...
                             'description',{'Localizations',    'Image size',   'Zoom',         'Number of time blocks'},...
                             'type',       {'array',            'array',        'array',        'array'},...
                             'dim_check',  {{[],[-1 2],[-1 3]}, {[],0,[1 2]},         {[],0},         {[],0},},...
                             'range_check',{'R',                'N+',           [eps Inf],      'N+'},...
                             'required',   {0,                  0,              0,              0},...
                             'default',    {'[]',               '[]',           1,             50}...
                              ),...
           'outparams',struct('name',{'frc_out'},...
                              'description',{'FRC curve'},...
                              'type',{'array'}...
                              )...
           );       
       
if nargin == 1
   s = varargin{1};
   if ischar(s) && strcmp(s,'DIP_GetParamList')
      frc_out = struct('menu','none');
      return
   end
end

try
   [positions, sz, zm, nblocks] = getparams(d,varargin{:});
catch
   if ~isempty(paramerror)
      error(paramerror)
   else
      error('Parameter parsing was unsuccessful.')
   end
end

% Check that position data are provided as input
if isempty(positions)
    if isempty(sz)
        error('No output image size can be defined.')
    else
        frc_out = zeros(1,ceil(max(sz)/sqrt(2)));
        return
    end
end

if isempty(zm)
   zm = 1;
end

% Check that the size of the binned images is provided
if isempty(sz)
   % Use a bounding box as output size
   sz = 1+round(zm*max(positions(:,1:2),[],1));
end

if isscalar(sz)
    sz = sz*[1 1];
end

% Check that a correct number of time blocks is specified
if isempty(nblocks)
    nblocks = 50;
end

if nblocks<2
    nblocks = 2;
    warning('postofrc:toofewblocks','Time series is split into 2 blocks.') 
end

% Assign localizations to half data sets
if size(positions,2) == 2
    if nblocks > size(positions,1)
        nblocks = size(positions,1);
        warning('postofrc:toomanyblocks','Number of time blocks is larger than the number of localizations.')
    end
    blocksel = randperm(nblocks,ceil(nblocks/2));                       % Block numbers in the first half data set              
    s = ceil((1:size(positions,1))/size(positions,1)*nblocks);          % Block numbers of the localizations
    s = ismember(s,blocksel);                                           % Localizations in the first half data set
else   
    % Sort blocks 
    maxt = max(positions(:,3));
    positions(:,3) = ceil(positions(:,3)/maxt*nblocks);                 % Go from time stamps to block numbers of the localizations                                      
    blocksel = randperm(nblocks,ceil(nblocks/2));                       % Block numbers in the first half data set
    s = ismember(positions(:,3),blocksel);                              % Localizations in the first half data set
end

% Bin images into 2 half data images
in1 = binlocalizations(positions(s,1:2),sz(1),sz(2),zm);
in2 = binlocalizations(positions(~s,1:2),sz(1),sz(2),zm);

% Compute FRC curve
frc_out = frc(in1,in2);
frc_out = double(frc_out);