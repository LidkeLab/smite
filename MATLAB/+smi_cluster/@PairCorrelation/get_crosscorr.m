function [C, r, c, dc, mask, rmax] = get_crosscorr(I1, I2, mask, rmax, flag)
% function [G, r, g, dg, mask] = get_crosscorr(I1, I2, mask, rmax, flag)
% calculates autocorrelation function for two dimensional images
%
% INPUTS
% I1 = First image to be autocorrelated
% I2 = Second image to be autocorrelated
% mask = region of interest, If none or [] specificed, user will be asked
%   to define.  To autocorreate the entire images, use mask = ones(size(I1))
% rmax = maximum r value to correlate in units of pixels. default is 100;
% flag = display flag.  insert 1 to display errorbar(r, g, dg) after
%   computation.
%
% OUTPUTS
% C = two dimensional cross-correlation function.  x and y values range between
%    -rmax:rmax
% r = radius values
% c = angularly averaged autocorrelation function.
% dc = errors on angularly averaged g
% mask = masked used for calculation
%
%
% Last updated 01.24.11 by Sarah Veatch.



if nargin<5, flag = 0; end  % flag for display
if (nargin<4 || isempty(rmax)), rmax=100; end  % distance of maximum correlation recorded
if (nargin<3 || isempty(mask)),    %% draw a mask if needed
    hIm = imshow(I1); axis equal tight off;
    %mask = roipoly;
    h = imrect;
    mask = createMask(h,hIm);
end

if sum(size(I1)==size(I2))<2,
    disp('images are not the same size')
    return
end

%mask = ones(size(I1)); end

N1 = sum(sum(I1.*mask));  % Total intensity within mask in I1
N2 = sum(sum(I2.*mask));  % Total intensity within mask in I2
A = sum(sum(mask));      % area of mask

I1 = double(I1);         % convert to double
I2 = double(I2);         % convert to double

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
C1 = A^2/N1/N2*real(fftshift(ifft2(fft2(I1.*mask,L1, L2).*conj(fft2(I2.*mask, L1, L2)))))./NP;
%G1 = A^2/N^2*real(fftshift(ifft2(abs(fft2(I1.*mask,L1, L2)).^2)))./NP; % 2D G(r) with proper normalization
C = imcrop(C1, [floor(L2/2+1)-rmax, floor(L1/2+1)-rmax, 2*rmax, 2*rmax]);  %only return valid part of G


xvals = ones(1, 2*rmax+1)'*(-rmax:rmax);    %map to x positions with center x=0
yvals = (-rmax:rmax)'*ones(1, 2*rmax+1);    %map to y positions with center y=0
zvals = C;

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
        c(j) = sum(m.*vv)/n2;               % the average G values in this bin
        dc(j) = sqrt(sum(m.*(vv-c(j)).^2))/n2; % the variance of the mean
    end
end

r = 0:rmax;

%end


if flag,
    r = 0:rmax;
    figure;errorbar(r, c, dc);
    axis tight
end

end
