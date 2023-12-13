%IMRES_LOCS   Compute resolution from a list of localizations
%
% The resolution is computed by taking blocks of localizations and assigning
% them randomly to half data sets. The localizations in these sets are
% binned into 2 images from which the FRC curve and subsequently the
% resolution are computed.
%
% SYNOPSIS:
%   [resolution_value frc_curve im_out resolution_high resolution_low resolution_t] = imres_locs(positions,size,zoomfactor,blocks,timefractions,reps,SR_pixelsize,positions_order, positions_units, show_im, show_frc,show_timefractions)
%
%   size
%      Size of the binned images
%   zoom
%      Pixels per unit distance in the positions list in the binned images
%   blocks
%      Number of time blocks into which the dataset is split.
%   timefractions
%      scalar: The acquisition time is split into timefractions blocks and
%              the resolution is computed for 1,2,...,timefractions blocks.
%      vector: The resolution is computed for each fraction of the
%              acquisition time specified in the vector.
%   reps
%      Number of times the FRC curve is computed for averaging
%   SR_pixelsize
%      Pixel size of the output image (in nm)
%   positions_order
%      Order of quantities in positions: 'xyt','xy', or 'txy'
%   positions_units
%      Units of quantities in positions: 'CCD pixels','SR pixels','nm', or 'um'
%   show_im
%      Display image of binned localizations?
%   show_frc
%      Display FRC curve?
%   show_timefractions
%      Display resolution vs time?
%
% DEFAULTS:
%   size : image edge corresponds to maximum y-coordinate
%   zoomfactor  = 1
%   blocks = 50
%   timefractions = 1
%   reps = 1 (no averaging)
%   SR_pixelsize = 10
%   positions_order = 'xyt'
%   positions_units = 'CCD pixels'
%   show_im = 1
%   show_frc = 1
%   show_timefractions = 1
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
% Robert Nieuwenhuizen & Bernd Rieger, Oct 2012


function varargout = imres_locs(varargin)

d = struct('menu','FRC resolution',...
    'display','Resolution from localizations',...
    'inparams',struct('name',       {'positions','sz', 'zoomfactor',   'blocks',                   'timefractions',                 'reps',                 'SR_pixelsize',                     'positions_order',  'positions_units',                      'show_im',          'show_frc',             'show_timefractions'},...
    'description',{'Localizations',    'Output image size',   'Zoom',  'Number of time blocks',    'Fractions of acquistion time',  'Number of averages',   'Superresolution pixel size (nm)',  'Input order',      'Input format',                         'Display image',    'Display FRC curve',    'Display resolution vs time'},...
    'type',       {'array',            'array',        'array',        'array',                    'array',                         'array',                'array',                            'option',           'option',                               'boolean',          'boolean',              'boolean'},...
    'dim_check',  {{[],[-1 2],[-1 3]}, {[],0,[1 2],[2 1]},         {[],0},         {[],0},         {[],0,[-1 1]},                   {[],0},                 0,                                  0,                  0,                                      0,                  0,                      0},...
    'range_check',{'R',                'N+',           [eps Inf],      'N+',                       'R+',                            'N+',                   [eps Inf],                          {'xy','xyt','txy'}, {'CCD pixels','SR pixels','nm','um'},   [],                 [],                     []},...
    'required',   {1,                  0,              0,              0,                          0,                               0,                      1,                                  0,                  0,                                      0,                  0,                      0},...
    'default',    {[],                 [],             1,              50,                         [],                              1,                      10,                                 'xyt',              'CCD pixels',                           1,                  1,                      1}...
    ),...
    'outparams',struct('name',{'resolution_value','frc_curve','im_out','resolution_high','resolution_low','resolution_t'},...
    'description',{'resolution value','FRC curve','Binned localizations','resolution - 1 std. dev.','resolution + 1 std. dev.','resolution vs time'},...
    'type',{'array','array','image','array','array','array'},...
    'suppress',{1,1,1,1,1,1}...
    )...
    );

if nargin == 1
    s = varargin{1};
    if ischar(s) & strcmp(s,'DIP_GetParamList')
        varargout{1} = d;
        return
    end
end

try
    [positions,sz,zoomfactor,nblocks,timefractions,reps,SR_pixelsize,positions_order, positions_units, show_im, show_frc,show_timefractions] = getparams(d,varargin{:});
catch
    if ~isempty(paramerror)
        error(paramerror)
    else
        error(firsterr)
    end
end

if isempty(zoomfactor)
    zoomfactor = 1;
end

%% Compute results

% Set right order of columns in positions
if strcmp(positions_order,'txy')
    positions = circshift(positions,[0 -1]);
end

% Convert units to superresolution pixels
switch positions_units
    case 'nm'
        positions(:,1:2) = positions(:,1:2)/SR_pixelsize;
        zoomfactor = 1;
    case 'um'
        positions(:,1:2) = positions(:,1:2)/1E3/SR_pixelsize;
        zoomfactor = 1;
end

% Compute results
[resolution_value, frc_curve, resolution_high, resolution_low, resolution_time] = postoresolution(positions, sz, zoomfactor,nblocks,timefractions,reps); % in super-resolution pixels
    
fprintf('Resolution value %2.1f +- %2.2f nm.\n', resolution_value*SR_pixelsize, (resolution_low-resolution_high)/2*SR_pixelsize);
fprintf('Resolution value %2.1f +- %2.2f superresolution pixels.\n', resolution_value, (resolution_low-resolution_high)/2);

%% Generate requested visualizations

% Output image
if isempty(sz)
    sz = 1+round(zoomfactor*max(positions(:,1:2),[],1));
end
if numel(sz) == 1
    sz = [sz sz];
end
im_out = mat2im(binlocalizations(positions, sz(1), sz(2), zoomfactor));
    
if show_im
    h=dipshow(gaussf(im_out,1));
    dipmapping(h,'colormap',hot)
end

% Plot FRC curve
if show_frc
    qmax = 0.5/(SR_pixelsize);
    
    figure
    hold on
    plot([0 qmax],[0 0],'k')
    plot(linspace(0,qmax*sqrt(2), length(frc_curve)), frc_curve,'-')
    plot([0 qmax],[1/7 1/7],'m')
    plot(1/(resolution_value*SR_pixelsize),1/7,'rx')
    plot(1/(resolution_value*SR_pixelsize)*[1 1],[-0.2 1/7],'r')
    hold off
    xlim([0,qmax]);
    ylim([-0.2 1.2])
    xlabel('Spatial frequency (nm^{-1})');
    ylabel('FRC')
end

% Show resolution value vs time

if isempty(timefractions) || all(timefractions==1)
    show_timefractions = 0;
end

if show_timefractions
    if isscalar(timefractions)
        timefractions = (1:timefractions)/timefractions;
    end
    
    figure
    plot(timefractions,resolution_time*SR_pixelsize)
    xlabel('Relative acquisition time');
    ylabel('Resolution (nm)');
    xlim([0 1])
end

%% Outputs
varargout{1} = resolution_value*SR_pixelsize;
varargout{2} = frc_curve;
varargout{3} = im_out;
varargout{4} = resolution_high*SR_pixelsize;
varargout{5} = resolution_low*SR_pixelsize;
varargout{6} = resolution_time*SR_pixelsize;
