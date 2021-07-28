function results = pair_correlation_Veatch(obj, SMD1, SMD2, correlation, fit)
% Pair correlation as originally written by Sarah L. Veatch.
%
% INPUTS:
%    obj             various properties used by the algorithms
%       BaseName        descriptive name for the results files
%       Fig_ext         figure extension
%       Font_props      [{'FontSize', 15, 'FontWeight', 'bold'}]
%       HistBinSize     pixel size (nm)
%       ResultsDir      directory to store results
%       ROI             [x_min, x_max, y_min, y_max]
%       Veatch_fit      ['exponential_and_gaussian'] fit model
%    SMD1 and SMD2   (x, y) coordinates of the two datasets (nm) in the format
%                    (1) SMD structures: SMD1.X, SMD1.Y, SMD2.X and SMD2.Y,
%                    (2) N x 2 array of coordinates.
%                    SMD2 is optional and if omitted or empty requests
%                    auto-correlation rather than cross-correlation.  Units are
%                    assumed to be nm.
%    correlation     'auto' or 'cross'
%    fit             [OPTIONAL] ['exponential_and_gaussian']
%                    'exponential_and_gaussian', 'exponential_and_cosine' or
%                    'exponential'
%
% OUTPUTS:
%    results         structure containing various results from the algorithm

   % Allow multiple input formats:
   if isstruct(SMD1) && isfield(SMD1, 'X') && isfield(SMD1, 'Y')
      XY1 = [SMD1.X, SMD1.Y];
   elseif ismatrix(SMD1) && size(SMD1, 2) == 2
      XY1 = SMD1;
   else
      error('SMD1 is not an SMD structure or an N x 2 matrix!');
   end

   if exist('SMD2', 'var') && ~isempty(SMD2)
      if isstruct(SMD2) && isfield(SMD2, 'X') && isfield(SMD2, 'Y')
         XY2 = [SMD2.X, SMD2.Y];
      elseif ismatrix(SMD2) && size(SMD2, 2) == 2
         XY2 = SMD2;
      else
         error('SMD2 is not an SMD structure or an N x 2 matrix!');
      end
   end

   base_name     = obj.BaseName;
   hist_bin_size = obj.HistBinSize;
   % Very inelegant, but I didn't want to break up this file into a bunch of
   % separate methods for a probably rarely if ever used functionality,
   % although it is nice to compare with pair_correlation results for
   % debugging purposes.
   property      = obj;

   if ~isempty(obj.ROI)
      x_min = obj.ROI(1);
      x_max = obj.ROI(2);
      y_min = obj.ROI(3);
      y_max = obj.ROI(4);
      ROI_size = min(obj.ROI([2, 4]) - obj.ROI([1, 3]));
   else
      if corr_type == 'C'
         x_min = min([XY1(:, 1); XY2(:, 1)]);
         x_max = max([XY1(:, 1); XY2(:, 1)]);
         y_min = min([XY1(:, 2); XY2(:, 2)]);
         y_max = max([XY1(:, 2); XY2(:, 2)]);
      else
         x_min = min(XY1(:, 1));
         x_max = max(XY1(:, 1));
         y_min = min(XY1(:, 2));
         y_max = max(XY1(:, 2));
      end
      ROI_size = min(x_max - x_min, y_max - y_min);
   end

   dataA(:, 1) = XY1(:, 1) - x_min;
   dataA(:, 2) = XY1(:, 2) - y_min;
   if exist('XY2', 'var') && ~isempty(XY2)
      dataB(:, 1) = XY2(:, 1) - x_min;
      dataB(:, 2) = XY2(:, 2) - y_min;
   else
      dataB = [];
   end

%% Filter the localisations
%rule = 1: precision or intensity
%rule = 2: both
% rule = 1;
% limits = [10 50]; %[min max]
% dataCol = 7; % variable to be used for filtering. Use a vector for two parameters eg [1 3]
% dataA = filter_localisations(dataA,dataCol,rule,limits);
% if ~isempty(dataB)
%     dataB = filter_localisations(dataB,dataCol,rule,limits);
% end
%% set some variables:
camPixSize = hist_bin_size; %pixel size on ccd in nm
originalX = round((x_max - x_min)/hist_bin_size);
originalY = round((y_max - y_min)/hist_bin_size); % image size in pixels
xScale = originalX .* camPixSize;
yScale = originalY .* camPixSize;
% nmPixSizeX = xScale / image_resolution;
% nmPixSizeY = yScale / image_resolution;
%nmPixSize = sqrt(nmPixSizeX^2 + nmPixSizeY^2); % pixel size in 2D histogram
nmPixSize = hist_bin_size;
image_resolution = [xScale/nmPixSize, yScale/nmPixSize]; % resolution for 2D histogram of localisation data
%image_resolution(1)
%% apply channel alignment?
transformation = []; %enter filename to apply transformation
calc_new = 0;
t_params = {transformation, calc_new};

%% set the type of correlation and the function to fit to the data
%%'auto' for auto-correlation, 'cross' for cross-correlation
if ~exist('correlation', 'var')
   correlation = 'auto';
   %correlation = 'cross';
end
if ~exist('fit', 'var')
   fit = obj.Veatch_fit;
   %fit = 'exponential_and_gaussian'; %name should match available fit functions
   %fit = 'exponential_and_cosine'; %name should match available fit functions
   %fit = 'exponential'; %name should match available fit functions
end

%radius = 1000; %in pixels
% Establish rmax as half the size of the ROI in pixels.
radius = round(ROI_size / (2 * hist_bin_size));

%% extract the (possibly filtered) x-y coordinates
% this step converts the coordinates pixels to nm; remove multiplication by
% camPixSize to work with data already in nm
Ax = dataA(:, 1);
Ay = dataA(:, 2);
Axy = [Ax Ay];
%Axy = [Ax Ay].*camPixSize;
if ~isempty(dataB)
    Bx = dataB(:, 1);
    By = dataB(:, 2);
    Bxy = [Bx By];
    %Bxy = [Bx By].*camPixSize;
else
    Bxy = [];
end

%% calculate the correlation and the fit
[correlation_data, params] = run_correlation_and_fit(Axy, Bxy, image_resolution, nmPixSize, t_params, [xScale yScale], correlation, fit, radius, base_name, property);

fprintf('Veatch %s-correlation using an %s fit for %s:', ...
        correlation, fit, base_name);
%correlation_data
params{1}
if numel(params) > 1
   params{2}
end

results.correlation_data = correlation_data;
results.params = params{1};

end

% -----------------------------------------------------------------------------

function [corrData, params] = run_correlation_and_fit(varargin) %(input1, input2, Xcol, Ycol, res, nm pix size, range, correlation type, fittype, radius)

    channel1 = varargin{1};
    channel2 = varargin{2};
    res = varargin{3};
    nmpixSize = varargin{4};
    t_params = varargin{5};
    range = varargin{6};
    Xrange = [0 range(1)];
    Yrange = [0 range(2)];
    correlation = varargin{7};
    fit = varargin{8};
    maxrad1 = varargin{9};
    base_name = varargin{10};
    property = varargin{11};

    %calculate histograms
    if ~isempty(channel2)
        density(:,:,1) = hist2d(channel1,res(1), res(2),Xrange,Yrange);
        density(:,:,2) = hist2d(channel2,res(1), res(2),Xrange,Yrange);
        density(:,:,3) = zeros(size(density(:,:,1),1),size(density(:,:,1),2));
    else
        density = hist2d(channel1,res(1), res(2),Xrange,Yrange);
    end
    N = {};
    data = {};
    N{1} = size(channel1,1);
    if isempty(channel2)
        numChannels = 1;
        data{1} = channel1;
        data{2} = [];
    else
        numChannels = 2;
        N{2} = size(channel2,1);
        if ~isempty(t_params{1}) && t_params{2} == 0
            tformData = load(t_params{1});
            TFORM = tformData.TFORM;
            data{1} = tformfwd(TFORM,channel1);
        elseif isempty(t_params{1}) && t_params{2} == 0
            data{1} = channel1;
        elseif isempty(t_params{1}) && t_params{2} == 1
            [in_points,base_points,~,TFORM] = transformChannels(nmpixSize,density);
            data{1} = tformfwd(TFORM,channel1);
        end
        data{2} = channel2;
    end
        %calculate histograms
    if numChannels == 2
        density(:,:,1) = hist2d(data{1},res(1), res(2),Xrange,Yrange);
        density(:,:,2) = hist2d(data{2},res(1), res(2),Xrange,Yrange);
        density(:,:,3) = zeros(size(density(:,:,1),1),size(density(:,:,1),2));
    else
        density = hist2d(data{1},res(1), res(2),Xrange,Yrange);
        cmax = 5;
        cmin = min(min(density));
    end

%   %display for ROI definition
%   LUT = RedMap;
%   r = LUT(:,2);
%   g = LUT(:,1);
%   b = LUT(:,3);
%   k = 1:numel(r);
%   map(k,:)=[r(k) g(k) b(k)];
%   figure;hIm = imagesc(Xrange,Yrange,density,[cmin,cmax]); axis equal tight off; colormap(map)
%   %mask = roipoly;
%   h = imrect;
%   pos = getPosition(h)
%   mask = createMask(h,hIm);
%   figure;imagesc(mask);axis equal tight off;colormap('gray')
%   figure;imagesc(density.*mask,[cmin,cmax]);axis equal tight off; colormap(map)
    params = {};
    corrData = repmat(struct('twoDcorr',[],'radius',[],'correlation',[],'error',[],'mask',[],'type',[],'L',[]),numChannels,1);
    L = {};
%     for ii = 1:numChannels
%         [corrData(ii).L,vq] = ripley(h,data{ii},maxrad1*nmpixSize);
%         fname = sprintf('clustermap0%d.tif',ii);
%         imwrite(uint16(vq),fname,'tif','compression','lzw')
%     end
    % the input parameters for the fit depend on the fit type
    switch fit
        case 'exponential_and_gaussian'
            P = [];
            P(1) = 800; % the decay of the exponential (in nm)
            P(2) = 10; %the amplitude of the exponential (y intercept)
            P(3) = 20; %sqrt(2)* the PSF of the image. (the half with of the gaussian) in nm
            P(4) = 1; % the surface density of the probe (in 1/um^2)
        case 'exponential_and_cosine'
            P = [];
            P(1) = 800; % the decay of the exponential (in nm)
            P(2) = 10; %the amplitude of the exponential (y intercept)
            P(3) = 20; %sqrt(2)* the PSF of the image. (the half with of the gaussian) in nm
        case 'exponential'
            P = [];
            P(1) = 800; % the decay of the exponential (in nm)
            P(2) = 10; %the amplitude of the exponential (y intercept)
    end

    %calculate correlation
    % Set to 1 to get a figure of g(r) with error bars.
    flag = 0;
    switch correlation
        case 'auto'
            for iChan = 1: numChannels
                Imsize1=min(size(density(:,:,1)));
                if Imsize1<1.25*maxrad1
                    maxrad1=round(Imsize1/1.25);
                end
                mask = ones(size(density(:,:,iChan)));
                [G, r, g, dg, maskout] = ...
                   smi_cluster.PairCorrelation.get_autocorr( ...
                      im2bw(density(:,:,iChan)), mask, maxrad1, flag);
                if isnan(g)
                    errordlg('auto-correlation calculation failed. try changing radius','modal');
                    return
                else
                    params{iChan} = fitData(r,g,dg,nmpixSize,fit,P,base_name,correlation,property);
                    corrData(iChan).twoDcorr = G;
                    corrData(iChan).radius = r;
                    corrData(iChan).correlation = g;
                    corrData(iChan).error = dg;
                    corrData(iChan).mask = maskout;
                    corrData(iChan).type = 'auto';
                end
            end
        case 'cross'
            Imsize1=min(size(density(:,:,1)));
            if Imsize1<1.25*maxrad1
                maxrad1=round(Imsize1/1.25);
            end

            if numChannels == 2
                mask = ones(size(density(:,:,1)));
                [C, r, c, dc, maskout] = ...
                   smi_cluster.PairCorrelation.get_crosscorr( ...
                      density(:,:,1), density(:,:,2), mask, maxrad1, flag);
            else
                errordlg('two channels are required to calculate cross-correlation','modal');
                return
            end

            if isnan(c)
                errordlg('cross-correlation calculation failed. try changing radius','modal');
                return
            else
                params{1} = fitData(r,c,dc,nmpixSize,fit,P,base_name,correlation,property);
                corrData(1).twoDcorr = C;
                corrData(1).radius = r;
                corrData(1).correlation = c;
                corrData(1).error = dc;
                corrData(1).mask = maskout;
                corrData(1).type = 'cross';
            end
    end
end

% -----------------------------------------------------------------------------

function params = fitData(x,y,err,pixSize,type,Pin,base_name,correlation,property)
%fit the data

    switch type
        case 'exponential_and_gaussian'
            params = repmat(struct('cluster_size',[],'magnitude',[],'density',[],'sigma',[]),1,1);
        case 'exponential_and_cosine'
            params = repmat(struct('cluster_size',[],'magnitude',[],'r0',[]),1,1);
        case 'exponential'
            params = repmat(struct('cluster_size',[],'magnitude',[]),1,1);
    end

    x = x .* pixSize;
    if ~isempty(property.Fig_ext)
       figure('Visible', 'off');
    else
       figure;
    end
    axes(property.Font_props{:});
    hold on
    errorbar(x(2:end), y(2:end), err(2:end), 'k.', 'LineWidth', 2)
    plot(x(2:end),ones(1, numel(x) - 1), 'b:', 'LineWidth', 3)
    switch type
        case 'exponential_and_gaussian'
            P0 = [Pin(1), Pin(2), Pin(3), Pin(4)];

            P = lsqcurvefit(@(P, r) exponential_and_gaussian(P, r), P0, x(2:end), y(2:end)-1);
            %P = lsqcurvefit('exponential_and_gaussian', P0, x(2:end), y(2:end)-1);

            %hold on
            plot(x, exponential_and_gaussian(P, x)+1, 'r', 'LineWidth', 3)
            %hold off

            %legend ('data', 'fit')

            params.cluster_size = P(1); %characteristic size of structure
            params.magnitude = P(2); % magnitude of clustering
            params.sigma = P(3)/sqrt(2); %in nm
            params.density = P(4); %in 1/um^2
        case 'exponential_and_cosine'
            P0 = [Pin(1), Pin(2), Pin(3)];
            P = lsqcurvefit(@(P, r) exponential_and_cosine(P, r), P0, x(2:end), y(2:end)-1);
            %P = lsqcurvefit('exponential_and_cosine', P0, x(2:end), y(2:end)-1);

            %hold on
            plot(x, exponential_and_cosine(P, x)+1, 'r', 'LineWidth', 3)
            %hold off

            %legend ('data', 'fit')

            params.cluster_size = P(1); %characteristic size of structure
            params.magnitude = P(2); % magnitude of clustering
            params.r0 = P(3); %in nm
        case 'exponential'
            P0 = [Pin(1), Pin(2)];
            P = lsqcurvefit(@(P, r) exponential(P, r), P0, x(2:end), y(2:end)-1);
            %P = lsqcurvefit('exponential', P0, x(2:end), y(2:end)-1);

            %hold on
            plot(x, exponential(P, x)+1, 'r', 'LineWidth', 3)
            %hold off

            %legend ('data', 'fit')

            params.cluster_size = P(1); %characteristic size of structure
            params.magnitude = P(2); % magnitude of clustering
    end
    legend(sprintf('%s-correlation data', correlation), 'g(r) Random', 'fit');
    axis tight
    xlabel('r (nm)');
    ylabel('g(r)');
    title(sprintf('%s (Veatch)', regexprep(base_name, '_', '\\_')));
    hold off
    if ~isdir(property.ResultsDir)
       mkdir(property.ResultsDir);
    end
    name = fullfile(property.ResultsDir, ...
                    sprintf('%s_%scorrV', base_name, correlation));
    if ~isempty(property.Fig_ext)
       print(['-d', property.Fig_ext], name);
       saveas(gcf, name);
    else
       saveas(gcf, name);
       delete(gcf);
    end
end

% =============================================================================

function [Hout Xbins Ybins] = hist2d(D, varargin) %Xn, Yn, Xrange, Yrange)
%HIST2D 2D histogram
%
% [H XBINS YBINS] = HIST2D(D, XN, YN, [XLO XHI], [YLO YHI])
% [H XBINS YBINS] = HIST2D(D, 'display' ...)
%
% HIST2D calculates a 2-dimensional histogram and returns the histogram
% array and (optionally) the bins used to calculate the histogram.
%
% Inputs:
%     D:         N x 2 real array containing N data points or N x 1 array
%                 of N complex values
%     XN:        number of bins in the x dimension (defaults to 20)
%     YN:        number of bins in the y dimension (defaults to 20)
%     [XLO XHI]: range for the bins in the x dimension (defaults to the
%                 minimum and maximum of the data points)
%     [YLO YHI]: range for the bins in the y dimension (defaults to the
%                 minimum and maximum of the data points)
%     'display': displays the 2D histogram as a surf plot in the current
%                 axes
%
% Outputs:
%     H:         2D histogram array (rows represent X, columns represent Y)
%     XBINS:     the X bin edges (see below)
%     YBINS:     the Y bin edges (see below)
%
% As with histc, h(i,j) is the number of data points (dx,dy) where
% x(i) <= dx < x(i+1) and y(j) <= dx < y(j+1). The last x bin counts
% values where dx exactly equals the last x bin value, and the last y bin
% counts values where dy exactly equals the last y bin value.
%
% If D is a complex array, HIST2D splits the complex numbers into real (x)
% and imaginary (y) components.
%
% Created by Amanda Ng on 5 December 2008

% Modification history
%   25 March 2009 - fixed error when min and max of ranges are equal.
%   22 November 2009 - added display option; modified code to handle 1 bin

    % PROCESS INPUT D
    if nargin < 1 %check D is specified
        error 'Input D not specified'
    end

    Dcomplex = false;
    if ~isreal(D) %if D is complex ...
        if isvector(D) %if D is a vector, split into real and imaginary
            D=[real(D(:)) imag(D(:))];
        else %throw error
            error 'D must be either a complex vector or nx2 real array'
        end
        Dcomplex = true;
    end

    if (size(D,1)<size(D,2) && size(D,1)>1)
        D=D';
    end

    if size(D,2)~=2;
        error('The input data matrix must have 2 rows or 2 columns');
    end

    % PROCESS OTHER INPUTS
    var = varargin;

    % check if DISPLAY is specified
    index = find(strcmpi(var,'display'));
    if ~isempty(index)
        display = true;
        var(index) = [];
    else
        display = false;
    end

    % process number of bins
    Xn = 20; %default
    Xndefault = true;
    if numel(var)>=1 && ~isempty(var{1}) % Xn is specified
        if ~isscalar(var{1})
            error 'Xn must be scalar'
        elseif var{1}<1 %|| ~isinteger(var{1})
            error 'Xn must be an integer greater than or equal to 1'
        else
            Xn = var{1};
            Xndefault = false;
        end
    end

    Yn = 20; %default
    Yndefault = true;
    if numel(var)>=2 && ~isempty(var{2}) % Yn is specified
        if ~isscalar(var{2})
            error 'Yn must be scalar'
        elseif var{2}<1 %|| ~isinteger(var{2})
            error 'Xn must be an integer greater than or equal to 1'
        else
            Yn = var{2};
            Yndefault = false;
        end
    end

    % process ranges
    if numel(var) < 3 || isempty(var{3}) %if XRange not specified
        Xrange=[min(D(:,1)),max(D(:,1))]; %default
    else
        if nnz(size(var{3})==[1 2]) ~= 2 %check is 1x2 array
            error 'XRange must be 1x2 array'
        end
        Xrange = var{3};
    end
    if Xrange(1)==Xrange(2) %handle case where XLO==XHI
        if Xndefault
            Xn = 1;
        else
            Xrange(1) = Xrange(1) - floor(Xn/2);
            Xrange(2) = Xrange(2) + floor((Xn-1)/2);
        end
    end

    if numel(var) < 4 || isempty(var{4}) %if XRange not specified
        Yrange=[min(D(:,2)),max(D(:,2))]; %default
    else
        if nnz(size(var{4})==[1 2]) ~= 2 %check is 1x2 array
            error 'YRange must be 1x2 array'
        end
        Yrange = var{4};
    end
    if Yrange(1)==Yrange(2) %handle case where YLO==YHI
        if Yndefault
            Yn = 1;
        else
            Yrange(1) = Yrange(1) - floor(Yn/2);
            Yrange(2) = Yrange(2) + floor((Yn-1)/2);
        end
    end

    % SET UP BINS
    Xlo = Xrange(1) ; Xhi = Xrange(2) ;
    Ylo = Yrange(1) ; Yhi = Yrange(2) ;
    if Xn == 1
        XnIs1 = true;
        Xbins = [Xlo Inf];
        Xn = 2;
    else
        XnIs1 = false;
        Xbins = linspace(Xlo,Xhi,Xn) ;
    end
    if Yn == 1
        YnIs1 = true;
        Ybins = [Ylo Inf];
        Yn = 2;
    else
        YnIs1 = false;
        Ybins = linspace(Ylo,Yhi,Yn) ;
    end

    Z = linspace(1, Xn+(1-1/(Yn+1)), Xn*Yn);

    % split data
    Dx = floor((D(:,1)-Xlo)/(Xhi-Xlo)*(Xn-1))+1;
    Dy = floor((D(:,2)-Ylo)/(Yhi-Ylo)*(Yn-1))+1;
    Dz = Dx + Dy/(Yn) ;

    % calculate histogram
    h = reshape(histc(Dz, Z), Yn, Xn);

    if nargout >=1
        Hout = h;
    end

    if XnIs1
        Xn = 1;
        Xbins = Xbins(1);
        h = sum(h,1);
    end
    if YnIs1
        Yn = 1;
        Ybins = Ybins(1);
        h = sum(h,2);
    end

    % DISPLAY IF REQUESTED
    if ~display
        return
    end

    [x y] = meshgrid(Xbins,Ybins);
    dispH = h;

    % handle cases when Xn or Yn
    if Xn==1
        dispH = padarray(dispH,[1 0], 'pre');
        x = [x x];
        y = [y y];
    end
    if Yn==1
        dispH = padarray(dispH, [0 1], 'pre');
        x = [x;x];
        y = [y;y];
    end

    surf(x,y,dispH);
    colormap(jet);
    if Dcomplex
        xlabel real;
        ylabel imaginary;
    else
        xlabel x;
        ylabel y;
    end

end

% -----------------------------------------------------------------------------

function y = exponential_and_gaussian(P, r)
xi = P(1);
A = P(2);
s = P(3);
rho = P(4)*1e-6;
if length(P)<5, C = 0; else C = P(5); end

A2 = 1/2/pi/s^2/rho;
y = A2*exp(-r.^2/2/s^2) + A*exp(-r/xi)+C;

end

% -----------------------------------------------------------------------------

function y = exponential_and_cosine(P, r)
xi = P(1);
A = P(2);
r0 = P(3);
if length(P)<4, C = 1; else C = P(4); end

y = C + A*exp(-r/xi).*cos((pi*r)/2/r0);

end

% -----------------------------------------------------------------------------

function y = exponential(P, r)
xi = P(1);
A = P(2);
if length(P)<3, C = 0; else C = P(3); end

y = A*exp(-r/xi)+C;

end
