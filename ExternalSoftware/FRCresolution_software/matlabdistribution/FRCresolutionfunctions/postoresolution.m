%POSTORESOLUTION   Compute resolution from a list of localizations
%
% The resolution is computed by taking blocks of localizations and assigning
% them randomly to half data sets. The localizations in these sets are
% binned into 2 images from which the FRC curve and subsequently the
% resolution are computed.
%
% SYNOPSIS:
%   [resolution frc_out resolution_high resolution_low resolution_t] = postoresolution(positions,size,zoomfactor,blocks,timefractions,reps)
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
%
% DEFAULTS:
%   size : image edge corresponds to maximum y-coordinate
%   zoomfactor  = 1
%   blocks = 50
%   timefractions = 1
%   reps = 1 (no averaging)
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

function varargout = postoresolution(varargin)

d = struct('menu','FRC resolution',...
    'display','Resolution from localizations',...
    'inparams',struct('name',       {'positions',        'sz', 'zm',    'blocks',                  'timefractions',                 'reps'},...
    'description',{'Localizations',    'Image size',   'Zoom',         'Number of time blocks',    'Fractions of acquistion time',  'Number of averages'},...
    'type',       {'array',            'array',        'array',        'array',                    'array',                         'array'},...
    'dim_check',  {{[],[-1 2],[-1 3]}, {[],0,[1 2]},         {[],0},         {[],0},                     {[],0,[-1 1]},                   0},...
    'range_check',{'R',                'N+',           [eps Inf],      'N+',                       'R+',                            'N+'},...
    'required',   {0,                  0,              0,              0,                          0,                               0},...
    'default',    {'[]',               '[]',           1,              50,                         '[]',                            1}...
    ),...
    'outparams',struct('name',{'resolution','frc_out','resolution_high','resolution_low','resolution_t'},...
    'description',{'Resolution','FRC curve','Upper bound','Lower bound','Resolution vs time'},...
    'type',{'array','array','array','array','array'}...
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
    [positions,sz,zm,nblocks,timefractions,reps] = getparams(d,varargin{:});
catch
    if ~isempty(paramerror)
        error(paramerror)
    else
        error(firsterr)
    end
end

if isempty(zm)
    zm = 1;
end

% Check that the size of the binned images is provided
if isempty(sz)
   % Use a square bounding box as output size
   sz = 1+round(zm*max(positions(:,1:2),[],1));
end

if reps > 1
%     matlabpool(4); %using the parallelization toolbox for faster computation
    frc_out = 0;
    res_tmp = zeros(reps,1);
    
%     parfor ii=1:reps
    for ii =1:reps
        %fprintf('  iteration %d/%d\n',ii,reps);
        % Calculate the FRC curve
        frc_out_tmp = postofrc(positions,sz,zm,nblocks);
        % Calculate the resolution from the FRC curve
        [res_tmp(ii)] = frctoresolution(frc_out_tmp,max(sz));
        frc_out = frc_out + frc_out_tmp;
    end
%     matlabpool close
    resolution = mean(res_tmp(res_tmp>0));
    frc_out  = frc_out./reps;
    % we do not divide by sqrt(reps) here as the repeats are correlated,
    % we keep it on the safe side this way (larger bounds).
    resolution_high = resolution-std(res_tmp);
    resolution_low = resolution+std(res_tmp);
else
    % Calculate the FRC curve
    frc_out = postofrc(positions,sz,zm,nblocks);
    % Calculate the resolution from the FRC curve
    [resolution resolution_high resolution_low] = frctoresolution(frc_out,max(sz));
end


% --- outputs ---
varargout{1} = resolution;
varargout{2} = frc_out;
if nargout > 2
    varargout{2} = frc_out;
    varargout{3} = resolution_high;
end
if nargout > 3
    varargout{4} = resolution_low;
end

% Check if resolution vs time needs to be calculated
if nargout < 5
    if ~isempty(timefractions);fprintf('timefractions given, but not requested as output.\n');end
    return;
end


%% Check that the variable timefractions is not trivial
if isempty(timefractions)
    varargout{5} = resolution;
    return
end
if (timefractions == 1)
    varargout{5} = resolution;
    return
end

% Check that positions is not empty
if isempty(positions)
    fprintf(' -- Could not find the resolution --\n')
    varargout{5} = (-1)*ones(length(timefractions,1));
    return
end

% Prepare calculation of resolution vs time
if isscalar(timefractions)
    if timefractions>1
        timefractions = (1:timefractions)/timefractions;
    end
else
    timefractions(timefractions>1) = 1;
end

if size(positions,2)==2
    N = length(positions); 
    T=0; %this is needed for parfor (altough never used)
else %3
    T = max(positions(:,3));
    N = 0; %this is needed for parfor (altough never used)
end

resolution_t = zeros(length(timefractions),1);
% Calculate resolution vs time
fprintf(' Computing resolution versus time ...\n')

% Check if the distributed computing toolbox is available
try 
    toolboxdir('distcomp');
    TB_distcomp=1;
catch
    TB_distcomp=0;
end

% Use for loop or parfor loop
if TB_distcomp

%     for nn = 1:length(timefractions)
    parfor nn = 1:length(timefractions)    
        % Select positions up to given fractions of time
        if size(positions,2)==2
            positions_tmp = positions(1:ceil(timefractions(nn)*N),:);
        else
            positions_tmp = positions(positions(:,3)<timefractions(nn)*T,:);
        end 
        if isempty(positions_tmp)
            resolution_t(nn) = -1;
            continue
        end

        res_tmp = zeros(reps,1);
        for ii=1:reps
            frc_out = postofrc(positions_tmp,sz,zm,nblocks);
            res_tmp(ii) = frctoresolution(frc_out,max(sz));
        end
        resolution_t(nn) = mean(res_tmp(res_tmp>0));
        %fprintf('res %d, %f\n',nn,resolution_t(nn));
    end
    % matlabpool close

else

    for nn = 1:length(timefractions)
        % Select positions up to given fractions of time
        if size(positions,2)==2
            positions_tmp = positions(1:ceil(timefractions(nn)*N),:);
        else
            positions_tmp = positions(positions(:,3)<timefractions(nn)*T,:);
        end 
        if isempty(positions_tmp)
            resolution_t(nn) = -1;
            continue
        end

        res_tmp = zeros(reps,1);
        for ii=1:reps
            frc_out = postofrc(positions_tmp,sz,zm,nblocks);
            res_tmp(ii) = frctoresolution(frc_out,max(sz));
        end
        resolution_t(nn) = mean(res_tmp(res_tmp>0));
        %fprintf('res %d, %f\n',nn,resolution_t(nn));
    end
    
end

varargout{5} = resolution_t;




