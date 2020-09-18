%QCORRECTION_IMS   Correct resolution from two images for spurious correlations
%
% SYNOPSIS:
%   [resolution_corr, resolution_uncorr, Q, frc_curve_corr frc_curve] = qcorrection_ims(in1, in2, meansig, stdsig, pixelsize, floorcor,show_frc)
%
%   pixelsize
%      Pixel size of the images (in nm)
%   show_frc
%      Display FRC curve?

%
% PARAMETERS:
%   meansig
%       Mean localization uncertainty (in SR pixels)
%   stdsig
%       St. dev. of localization uncertainties (in SR pixels)
%   pixelsize  
%      Pixel size of the images (in nm)
%   floorcor
%      Take into account noise floor on the numerator of the FRC?
%
% OUTPUT:
%   resolution_corr
%      Resolution value after correction for spurious correlations
%   resolution_uncorr
%      Resolution value without correction
%   Q                   
%      Estimate for the number of times an emitter is localized on average 
%      assuming Poisson statistics for the localizations per emitter
%  frc_curve_corr
%      FRC curve corrected for spurious correlations
%  frc_curve            
%      FRC curve without correction
%
% NOTES:
%   Non-square images are zero padded to make them square.
%
% SEE ALSO:
%  binlocalizations, frc, frctoresolution
%
% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen & Bernd Rieger, Dec 2012

function varargout = qcorrection_ims(varargin)

d = struct('menu','FRC resolution',...
           'display','Q-corrected resolution from images',...
           'inparams',struct('name',       {'in1',      'in2',          'meansig',          'stdsig',                   'SR_pixelsize',     'floorcorr',                'show_frc'},...
                             'description',{'Image 1',  'Image 2',      'Loc. unc.',        'Std. dev. of loc. unc.'    'Pixel size (nm)',  'Correct for noise floor',  'Display figures'},...
                             'type',       {'image',    'image',        'array',            'array',                    'array',            'boolean',                  'boolean'},...
                             'dim_check',  {2,          2,              0,                  0,                          0,                  0,                          0},...
                             'range_check',{[],         [],             [eps Inf],          'R+',                       [eps Inf],          [],                         []},...
                             'required',   {1,          1,              1,                  1,                          1,                  0,                          0},...
                             'default',    {'in1',      'in2',          1,                  0,                          10,                 0,                          0}...
                              ),...
           'outparams',struct('name',{'resolution_corr','resolution_value','Q','frc_curve_corr','frc_curve'},...
                              'description',{'Resolution (corrected)','Resolution (uncorrected)','Q','FRC curve (corrected)','FRC curve (uncorrected)'},...
                              'type',{'array','array','array','array','array'},...
                              'suppress',{1,1,1,1,1}...
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
   [in1, in2, meansig, stdsig, SR_pixelsize, floorcorr, show_frc] = getparams(d,varargin{:});
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

%% Make images square
% Compute mask in x-direction
nfac = 8;                                                   % Image width / Width of edge region
x_im = xx(sz(1),sz(2))/sz(1);
mask = 0.5-0.5*cos(pi*nfac*x_im);          
mask(abs(x_im)<((nfac-2)/(nfac*2))) = 1;

% Check that input images are square and mask
if sz(1) == sz(2)
    mask = mask*rot90(mask);
    
    % Mask input images
    in1 = mask*in1;
    in2 = mask*in2;
    
else
    warning('qcorrection_ims:nonsquare','Images are not square.');
   
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
    in1 = extend(in1,[max(sz) max(sz)]);
    in2 = extend(in2,[max(sz) max(sz)]);
end

%% Prepare calculation of Q
Qnorm = (1/mean(in1)+1/mean(in2));

%% Fourier transform input images
in1 = ft(in1);
in2 = ft(in2);

%% Calculate FRC numerator
frcnum = double(real(radialmean(in1.*conj(in2))));

sz = size(in1,1);
q = (0:(length(frcnum)-1))./sz;

% find logmin and minloc for range of scaled sigs in order to be able to
% find average noise floor in that region, the noisefloor varies steeply in
% the plateau-region depending on where the minimum is on the plateau (to
% the left or to the right) the noisefloor-correction will be too large or
% too small, respectively.

allscales = [0.8 1.2 1.0];
allminloc = zeros(3,1);
alllogmin = zeros(3,1);

for jsc=1:numel(allscales)
scale = allscales(jsc);

% Calculate exponential decay function for Gaussian distribution of sigs
stdfac = 1+8*pi^2*(scale*stdsig)^2*(q.^2);
exp_decay = exp(-4*pi^2*(scale*meansig)^2*(q.^2)./stdfac)./sqrt(stdfac);

% Find smoothing kernel size
tmplog = real(log(frcnum) - log(exp_decay) -2*log(sinc(q)));
tmplog = tmplog(~isnan(tmplog));
tmplog = tmplog(~isinf(tmplog));
d = find(tmplog>(1+tmplog(1)),1);

if ((isempty(d) || d>sz/2))
    d = sz/2;
end

% Smooth FRC numerator
tmplog = real(log(frcnum) - log(exp_decay) -2*log(sinc(q)));
smoothlog = cfsmooth(real(log(frcnum) - log(exp_decay) -2*log(sinc(q))),d/5,'rloess');
if any(smoothlog==0)
    smoothlog(smoothlog==0) = tmplog(smoothlog==0);
end

% Find minimum of logarithm
% smoothlog(1:99) = 1E10;
% smoothlog(801:end) = 1E10;
[logmin minloc] = min(smoothlog(1:ceil(sz/2)));
allminloc(jsc) = minloc;
alllogmin(jsc) = logmin;

end

% Correct for contribution from the noise floor
nr = double(radialsum(ones(sz)));
noisefac = mean(log(frcnum(round(sz/2):end).^2.*nr(round(sz/2):end)));
noiselog = real(0.5*log(exp(noisefac)./nr)) - log(exp_decay) -2*log(sinc(q));

meannoiselog = (noiselog(allminloc(1))+noiselog(allminloc(2)))/2;
if floorcorr
  if logmin>meannoiselog
    logmin = log(exp(logmin)-exp(meannoiselog));
  else
    logmin = -Inf;
  end
end

% Correct for spurious correlations
frcnum_corr = exp(logmin).*exp_decay.*(sinc(q)).^2;

% Estimate number of localizations per emitter
Q = exp(logmin)*Qnorm;
fprintf('Estimated Q: %2.2f.\n',Q);

% Corrected frc
frcdenom = double(sqrt(real(radialmean(abs(in1).^2).*radialmean(abs(in2).^2))));
frccurve = frcnum./frcdenom;
frccurve_corr = (frcnum-frcnum_corr)./(frcdenom +frcnum_corr);


%% Plot results
if show_frc
    figure
    hold on
    plot(q/SR_pixelsize,frccurve);
    plot(q/SR_pixelsize,frccurve_corr,'r')
    plot([0 1/SR_pixelsize],[0 0],'k');
    xlim([0 0.5/SR_pixelsize])
    xlabel('Spatial frequency (nm^{-1})')
    ylabel('FRC')
    legend('Normal','Corrected for spurious correlations')
    hold off

    figure
    hold on
    % semilogy(q/ps, exp(log(Qnorm)+real(log(frcnum)) - log(exp_decay) -2*log(sinc(q))),'k',...
    %          q/ps, (log(Qnorm)+smoothlog)/log(10),'r','Linewidth',2,...
    %          q/ps, (log(Qnorm)+noiselog)/log(10),'g','Linewidth',2,...
    %          q/ps, (log(Qnorm)+smoothlog(minloc)*ones(size(q)))/log(10),'m','Linewidth',2)
    plot(q/SR_pixelsize, (log(Qnorm)+real(log(frcnum)) - log(exp_decay) -2*log(sinc(q)))/log(10),'k')
    plot(q/SR_pixelsize, (log(Qnorm)+smoothlog)/log(10),'r')
    % plot(q/ps, (log(Qnorm)+noiselog)/log(10),'g','Linewidth',2);
    plot(q/SR_pixelsize, (log(Qnorm)+smoothlog(minloc)*ones(size(q)))/log(10),'m')
    % plot(q(allminloc(:))/ps, (log(Qnorm)+smoothlog(allminloc(:)))/log(10),'ro');
    xlim([0 0.5/SR_pixelsize])
    ylim([min(log(frcnum(frcnum>0))) log(log(Qnorm)+frcnum(1))])
    ylim([-2 3])
    % % hleg = legend('normal','smooth','noise floor','plateau','minimum','Position','SouthEast');
    legend('Raw data','Smoothed data','Found plateau');
    hold off
    xlabel('Spatial frequency (nm^{-1})')
    ylabel('^{10}log(scaled FRC numerator)')
    set(gca,'YTick',[-2 -1 0 1 2 3])
    box on

    figure
    hold on
    plot(q/SR_pixelsize, real(log(frcnum)))
    plot(q/SR_pixelsize, real(log(frcdenom)),'k')
    plot(q/SR_pixelsize, real(log(frcnum_corr)),'r')
    plot(q/SR_pixelsize, real(log(frcnum-frcnum_corr)),'m')
    plot(q/SR_pixelsize, real(0.5*log(exp(noisefac)./nr)),'g');
    xlim([0 0.5/SR_pixelsize])
    ylim([min(log(frcnum(frcnum>0))) log(frcnum(1))])
    legend('FRC numerator','FRC denominator','Estimated spurious term','Corrected numerator','Noise floor')
    hold off
    xlabel('Spatial frequency (nm^{-1})')
    ylabel('log(FRC numerator)')
end

%% Calculate resolutions
resolution_uncorrr = frctoresolution(frcnum./frcdenom,sz);
resolution_corr = frctoresolution((frcnum-frcnum_corr)./(frcdenom +frcnum_corr),sz);

varargout{1} = resolution_corr;
varargout{2} = resolution_uncorrr;
varargout{3} = Q;
varargout{4} = frccurve_corr;
varargout{5} = frccurve;