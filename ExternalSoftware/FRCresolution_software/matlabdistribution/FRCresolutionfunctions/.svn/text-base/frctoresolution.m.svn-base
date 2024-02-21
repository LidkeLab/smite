%FRCTORESOLUTION   Compute the resolution from a FRC curve
%
% SYNOPSIS:
%   [resolution resolution_high resolution_low] = binlocalizations(frc_in,size)
%
% PARAMETERS:
%   size
%     Size in pixels of the half data images from which frc_in was obtained

% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen, Oct 2012

function varargout = frctoresolution(frc_in,sz)

% Check the format of frc_in
if ~isa(frc_in,'numeric')
    error('The FRC curve does not have the right data format.')
end
if ~isvector(frc_in)
    error('The FRC curve has the wrong dimensionality')
end
if isrow(frc_in)
    frc_in = frc_in';
end
    % Check the format of sz
if ~isa(sz,'numeric')
    error('The specified image size does not have the right data format.')
end
if ~isscalar(sz)
    if numel(sz)~= 2
        error('The specified image size does not have the right data format.');
    else
        sz = max(sz);
    end
end
if (sz ~= round(sz)) || ~(sz > 0) 
    error('The image size must be a positive integer.')
end

% Check that the curvefit toolbox function smooth exists
TB_curve=0;
try 
    TB_d=toolboxdir('curvefit');
    TB_curve=1;
catch
    warning('Curvefit toolbox not available. Using another not optimal smoothing method for FRC.')
end

%% Check for undesirable values in FRC curve
frc_in = real(frc_in);
% frc_in(isnan(frc_in)) = 0;
frc_in(frc_in>1) = 1;
frc_in(frc_in<-1) = -1;

%% Calculate the number of pixels in a Fourier ring
% Number of pixels in a Fourier ring
nr = double(radialsum(ones(sz)))';       

if length(nr) ~= length(frc_in)
   error('The specified image size could not have resulted in this FRC curve.'); 
end

% Calculate the effective width of pixels in a Fourier ring
if nargout >1
    n1 = find(frc_in<0, 1,'first' );
    frc_tmp = frc_in;
    frc_tmp(isnan(frc_tmp)) = 0;
    frc_tmp(isinf(frc_tmp)) = 0;
    ne_inv = mean(frc_tmp(n1:end).^2.*nr(n1:end));
end

%% Smoot the FRC curve 
% Least squares interpolation for curve smoothing
sspan = ceil(sz/20);      % Smoothing span
if (sz/20)<5
    sspan = 5;
end
sspan = sspan  + (1-mod(sspan,2));

if TB_curve
   p = pwd; % hack to avoid the function shadwoing by smooth from dip_image
   cd([TB_d filesep 'curvefit'])
   frc_in = double(smooth(frc_in,sspan,'loess'));        
   cd(p)
else
   frc_in = [double(gaussf(frc_in,.9))]';
end
q = (0:(length(frc_in)-1))'/sz;                   % Spatial frequencies

%% Calculate intersections between the FRC curve and the threshold curve
% isects = polyxpoly(q,frc_in,q,thresholdcurve);
thresholdcurve = 1/7*ones(size(frc_in));
isects = isect(q,frc_in,thresholdcurve);

%% Find first intersection to obtain the resolution

% Throw away intersections at frequencies beyond the Nyquist frequency
isects = isects(isects<0.5);

if ~isempty(isects)
    % Find the first intersection where the FRC curve is decreasing
    isect_inds = 1+floor(sz*isects);    % Indices of the intersections
    for ii = 1:length(isect_inds)
        isect_ind = isect_inds(ii);
        if frc_in(isect_ind+1) < frc_in(isect_ind)
            resolution = 1/isects(ii);
            break
        end
    end
end

if ~exist('resolution','var')
    resolution = -1;
end

if resolution == -1
    varargout{1} = NaN;
    varargout{2} = 0;
    varargout{3} = 0.5;
    fprintf(' -- Could not find the resolution --\n')
    return
end

if nargout < 2
    varargout{1} = resolution;
    return
end

%% Calculate the resolution uncertainty

% Variance of the unsmoothed FRC values
frc_var = (1+2*frc_in-frc_in.^2).*(1-frc_in).^2*ne_inv./nr;

% Indices of FRC values inside the smoothing window
indexvec = (ceil(isect_ind-sspan/2)):(floor(isect_ind+sspan/2));    
indexvec = indexvec(indexvec>0);
indexvec = indexvec(indexvec<=length(frc_in));

% Calculate the covariance matrix of the LMS fit parameters.
S = diag(frc_var(indexvec));    % Approximate covariance matrix of the FRC values
X = cat(2,ones(size(indexvec')),q(indexvec),q(indexvec).^2);
C = inv(X'*X);
C = C*X'*S*X*C;

% Estimate the variance of FRC(fres) resulting from the least mean
% squares fit. The factor ne_inv arises to correct for the
% correlations between neighbouring FRC values.
frc_newvar = ne_inv*(1./[1 resolution resolution^2])*C*(1./[1; resolution; resolution^2]);

if isnan(frc_newvar) || isinf(frc_newvar)
    frc_newvar = frc_var(isect_ind);
end

% Calculate the crossing between FRC + 1*stdev(FRC) and threshold
frc_high = frc_in+sqrt(frc_newvar);
isects_high = isect(q,frc_high,thresholdcurve);
isects_high = isects_high(isects_high<0.5);
% Calculate the crossing between FRC - 1*stdev(FRC) and threshold 
frc_low = frc_in-sqrt(frc_newvar);
isects_low = isect(q,frc_low,thresholdcurve);
isects_low = isects_low(isects_low<0.5);

% Find resolution + 1 standard deviation
if ~isempty(isects_high)
    isect_inds = 1+floor(sz*isects_high);
    for ii = 1:length(isect_inds)
        isect_ind = isect_inds(ii);
        if frc_high(isect_ind+1) < frc_high(isect_ind)
            resolution_high = 1/isects_high(ii);
            break
        end
    end
end
if ~exist('resolution_high','var')
    resolution_high = 0;
end

% Find resolution + - standard deviation
if ~isempty(isects_low)
    isect_inds = 1+floor(sz*isects_low);
    for ii = 1:length(isect_inds)
        isect_ind = isect_inds(ii);
        if frc_low(isect_ind+1) < frc_low(isect_ind)
            resolution_low = 1/isects_low(ii);
            break
        end
    end
end
if ~exist('resolution_low','var')
    resolution_low = 0.5;
end

%% Write outputs
varargout{1} = resolution;
varargout{2} = resolution_high;
varargout{3} = resolution_low;
