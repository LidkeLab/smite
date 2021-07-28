function [G, r, g, dg, mask, rmax] = get_autocorr(I1, mask, rmax, flag)
% function [G, r, g, dg, mask] = get_autocorr(I1, mask, rmax, flag)
% calculates autocorrelation function for two dimensional images
%
% INPUTS
% I1 = image to be autocorrelated
% mask = region of interest, If none or [] specificed, user will be asked
%   to define.  To autocorreate the entire images, use mask = ones(size(I1))
% rmax = maximum r value to correlate in units of pixels. default is 100;
% flag = display flag.  insert 1 to display errorbar(r, g, dg) after
%   computation.
%
% OUTPUTS
% G = two dimensional correlation function.  x and y values range between
%    -rmax:rmax
% r = radius values
% g = angularly averaged autocorrelation function.
% dg = errors on angularly averaged g
% mask = masked used for calculation
%
% NOTE: G(r=0) is just the dot product of the image.  For display purposes,
% G(r=0) is set to zero in the 2D autocorrelation output.  g(r=0) [g(1)]
% retains the proper value.
%
% Last updated 01.26.10 by Sarah Veatch.



if nargin<4, flag = 0; end  % flag for display
if (nargin<3 || isempty(rmax)), rmax=100; end  % distance of maximum correlation recorded
if (nargin<2 || isempty(mask)),    %% draw a mask if needed
    imagesc(I1); axis equal tight off;
    mask = roipoly;
end

%mask = ones(size(I1)); end

N = sum(sum(I1.*mask));  % number of particles within mask
A = sum(sum(mask));      % area of mask

I1 = double(I1);         % convert to double

L1 = size(I1, 1)+rmax; % size of fft2 (for zero padding)
L2 = size(I1, 2)+rmax; % size of fft2 (for zero padding)

% Adjust rmax if it is too big and would cause problems cropping C1:
%    Need 2*rmax + 1 >= min(size(I1, 1), size(I1, 2)) + rmax
% This fix is for rectangular regions which are treated somewhat differently in
% get_corr, so rectangular regions may produce different results in the two
% routines---use with caution!  (Square region results are identical.)
Lmin = min(L1, L2);
if 2*rmax + 1 > Lmin
   fprintf('rmax adjusted from %.1f to ', rmax);
   rmax = min(size(I1, 1), size(I1, 2)) - 1;
   fprintf('%.1f\n', rmax);

   L1 = size(I1, 1)+rmax; % size of fft2 (for zero padding)
   L2 = size(I1, 2)+rmax; % size of fft2 (for zero padding)
end

NP = real(fftshift(ifft2(abs(fft2(mask, L1, L2)).^2))); % Normalization for correct boundary conditions
G1 = A^2/N^2*real(fftshift(ifft2(abs(fft2(I1.*mask,L1, L2)).^2)))./NP; % 2D G(r) with proper normalization
G = imcrop(G1, [floor(L2/2+1)-rmax, floor(L1/2+1)-rmax, 2*rmax, 2*rmax]);  %only return valid part of G


xvals = ones(1, 2*rmax+1)'*(-rmax:rmax);    %map to x positions with center x=0
yvals = (-rmax:rmax)'*ones(1, 2*rmax+1);    %map to y positions with center y=0
zvals = G;

[theta,r,v] = cart2pol(xvals,yvals, zvals);  % convert x, y to polar coordinates

Ar = reshape(r,1, (2*rmax+1)^2);
Avals = reshape(v,1, (2*rmax+1)^2);
[rr,ind] = sort(Ar);                         % sort by r values
vv = Avals(ind);                             % reindex g
r = 0:floor(max(rr));                        % the radii you want to extract
[n bin] = histc(rr, r-.5);                   % bin by radius

for j = 1:rmax+1;                            % now get averages
    m = bin==j;
    n2 = sum(m);                             % the number of pixels in that bin
    if n2==0, vals(j)=0; er(j)=0;            % if no bins, no data
    else
        g(j) = sum(m.*vv)/n2;               % the average G values in this bin
        dg(j) = sqrt(sum(m.*(vv-g(j)).^2))/n2; % the variance of the mean
    end
end

r = 0:rmax;

%end

G(rmax+1, rmax+1) = 0;

if flag,
    r = 0:rmax;
    errorbar(r(2:length(r)), g(2:length(r)), dg(2:length(r)));
    axis tight
end

end
