%WKK   Sum values along lines perpendicular to each position vector [x y]
%
% A line is defined through each vector which is perpendicular to the line
% from that pixel to the center pixel of the image. The sum of the pixel 
% values along that line is taken as outfput.
%
% SYNOPSIS:
%   out = wkk(in)
% 
% NOTES:
%  The center pixel of in is defined as: floor((size(in)-1)/2) 

% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen, Oct 2012

function out = wkk(in)

in = real(in);

% Check that the input is a square 2D image.
if ~isa(in,'dip_image')
    error('Input must be an image.');
end

if ~(ndims(in) == 2)
    error('Input image must be 2D.')
end
    
% Radon transform image
t = 0:1:180;                    % Angles
[c,xp] = radon(im2mat(in),-t);

% Find xy-coordinates of radon transform values
[t,xp] = ndgrid(t./180*pi,xp);
[x,y] = pol2cart(t,xp);

% Compute grid of output coordinates
xc = (-ceil((size(in,1)-1)/2):floor((size(in,1)-1)/2));
yc = (-ceil((size(in,2)-1)/2):floor((size(in,2)-1)/2))';

% Check if a warning will be displayed when calling TriScatteredInterp
s = warning('query', 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId');
if strcmp(s.state,'on')
    warning('off', 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId');
end

% Interpolate the radon transform data on the gridpoints (xc,yc)
Int_obj = TriScatteredInterp(reshape(x,[numel(x) 1]),reshape(y,[numel(x) 1]),reshape(c',[numel(x) 1]));
out = mat2im(Int_obj(repmat(xc,[length(yc) 1]),repmat(yc,[1 length(xc)])));
out(isnan(out)) = 0;

% Return the state of the warning from griddata to its original state
if strcmp(s.state,'on')
    warning('on', 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId');
end
