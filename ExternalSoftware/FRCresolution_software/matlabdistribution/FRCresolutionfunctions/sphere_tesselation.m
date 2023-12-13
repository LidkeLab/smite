%SPHERE_TESSELATION   Generate 4^N sampling angles
%
% SYNOPSIS:
%   [polar_angle azimuthal_angle] = sphere_tesselation(N)
%
% NOTES:
%   This function is for generating a tesselation of a sphere by iteratively
%   dividing triangles in 4 equal triangles. Starting from 8 octants defined
%   by +/-e_x, +/-e_y, +/-e_z we get in n steps 2*4^n triangles. Taking the
%   center of gravity of the three corners gives 2*4^n directions. Rotating
%   over R to put one of these on +e_z and using those vectors that have
%   n_z>0 should give 4^n indpendent directions covering the half-sphere.

% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Sjoerd Stallinga, Oct 2012

function [polang,aziang] = sphere_tesselation(maxlevel)
% Avoid being in menu
if nargin==1
   if ischar(varargin{1}) & strcmp(varargin{1},'DIP_GetParamList')
      polang = struct('menu','none');
      return
   end
end
if ~isnumeric(maxlevel)
    error('Input must be a positive integer.');
end

if ~isscalar(maxlevel)
    error('Input must be a positive integer.');
end

if ~(maxlevel>0)
    error('Input must be a positive integer.');
end

alltriangles = cell(maxlevel,2*4^maxlevel);

% initial level
plx = [1,0,0];
mnx = [-1,0,0];
ply = [0,1,0];
mny = [0,-1,0];
plz = [0,0,1];
mnz = [0,0,-1];
triangleinit = cell(8,1);
triangleinit{1} = [plx;ply;plz];
triangleinit{2} = [ply;mnx;plz];
triangleinit{3} = [mnx;mny;plz];
triangleinit{4} = [mny;plx;plz];
triangleinit{5} = [plx;ply;mnz];
triangleinit{6} = [ply;mnx;mnz];
triangleinit{7} = [mnx;mny;mnz];
triangleinit{8} = [mny;plx;mnz];
for jj=1:8
  alltriangles{1,jj} = triangleinit{jj};
end

% loop over all levels
for ilevel = 2:maxlevel
  for jj=1:2*4^(ilevel-1)
    allvecs = alltriangles{ilevel-1,jj};
    a = allvecs(1,:);
    b = allvecs(2,:);
    c = allvecs(3,:);
    ab = a+b;
    norm = sqrt(sum(ab.^2));
    ab = ab/norm;
    bc = b+c;
    norm = sqrt(sum(bc.^2));
    bc = bc/norm;
    ca = c+a;
    norm = sqrt(sum(ca.^2));
    ca = ca/norm;
    alltriangles{ilevel,4*jj-3} = [a;ab;ca];
    alltriangles{ilevel,4*jj-2} = [b;bc;ab];
    alltriangles{ilevel,4*jj-1} = [c;ca;bc];
    alltriangles{ilevel,4*jj} = [ab;bc;ca];
  end
end

% find centers of gravity of all triangles
alldirections = zeros(2*4^maxlevel,3);
for jj=1:2*4^maxlevel
  allvecs = alltriangles{maxlevel,jj};
  a = allvecs(1,:);
  b = allvecs(2,:);
  c = allvecs(3,:);
  abc = a+b+c;
  norm = sqrt(sum(abc.^2));
  abc = abc/norm;
  alldirections(jj,:) = abc;
end

% generate rotation to get one direction along +e_z
cosphi = 1/sqrt(2.0);
sinphi = 1/sqrt(2.0);
costheta = 1/sqrt(3.0);
sintheta = sqrt(2.0/3.0);
Rz = [cosphi,-sinphi,0;sinphi,cosphi,0;0,0,1];
Rx = [1,0,0;0,costheta,-sintheta;0,sintheta,costheta];
R = Rx*Rz;
for jj=1:2*4^maxlevel
  alldirections(jj,:) = alldirections(jj,:)*R';
end

% restrict to vectors with n_z>0
epsy = 1e-6;
alldirections = alldirections(alldirections(:,3)>-epsy,:);

% transform to spherical coordinates
[aziang,polang,radius] = cart2sph(alldirections(:,1),alldirections(:,2),alldirections(:,3));
polang = pi/2-polang;

% % make plot
% figure;
% u = alldirections(:,1);
% v = alldirections(:,2);
% w = alldirections(:,3);
% allzeros = zeros(size(alldirections));
% x = allzeros(:,1);
% y = allzeros(:,2);
% z = allzeros(:,3);
% quiver3(x,y,z,u,v,w);

end
