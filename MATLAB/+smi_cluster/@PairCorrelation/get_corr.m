function [C, r, c, dc, rmax] = get_corr(n_ROIs, rmax, II1, II2)
% function [G, r, g, dg] = get_crosscorr(I1, I2, rmax)
% calculates autocorrelation function for two dimensional images
%
% INPUTS
% nROIs = number of ROIs to be correlated
% II1 = cell array of ROI images in first  image to be autocorrelated
% II2 = cell array of ROI images in second image to be autocorrelated
% rmax = maximum r value to correlate in units of pixels.
%
% OUTPUTS
% C = two dimensional cross-correlation function.  x and y values range between
%    -rmax:rmax
% r = radius values
% c = angularly averaged autocorrelation function.
% dc = errors on angularly averaged g
% rmax = adjusted input rmax value
%
% NOTE: G(r=0) is just the dot product of the image.  For display purposes,
% G(r=0) is set to zero in the 2D autocorrelation output.  g(r=0) [g(1)]
% retains the proper value.
%
% Last updated 01.24.11 by Sarah Veatch.

% Modified by
%    Michael J. Wester (2021)

if exist('II2', 'var')
   corr_type = 'C';
else
   corr_type = 'A';
end

L1min = 1e+10;
L2min = 1e+10;
L1max = -1;
L2max = -1;
for i = 1 : n_ROIs
   L1 = size(II1{i}, 1)+rmax; % size of fft2 (for zero padding)
   L2 = size(II1{i}, 2)+rmax; % size of fft2 (for zero padding)
   L1min = min(L1min, L1);
   L2min = min(L2min, L2);
   L1max = max(L1max, L1);
   L2max = max(L2max, L2);
end
Lmax_min = min(L1max, L2max);

% Adjust rmax if it is too big and would cause problems cropping C1:
%    Need floor(Lmax_min/2+1) - rmax >= 1
if floor(Lmax_min/2) < rmax
   fprintf('rmax adjusted from %.1f to ', rmax);
   rmax = floor(Lmax_min/2);
   fprintf('%.1f\n', rmax);
end

sum_A  = 0;
sum_N1 = 0;
sum_N2 = 0;
sum_FF = 0;
sum_NP = 0;
for i = 1 : n_ROIs
   I1 = II1{i};
   if corr_type == 'C'
      I2 = II2{i};
      if sum(size(I1)==size(I2))<2,
         disp('images are not the same size')
         return
      end
   end

   N1 = sum(sum(I1));        % Total intensity in I1
   if corr_type == 'C'
      N2 = sum(sum(I2));     % Total intensity in I2
   end
   A = prod(size(I1));       % area of mask
   mask = ones(size(I1));

   I1 = double(I1);          % convert to double
   if corr_type == 'C'
      I2 = double(I2);       % convert to double
   end

   %L1 = size(I1, 1)+rmax;    % size of fft2 (for zero padding)
   %L2 = size(I1, 2)+rmax;    % size of fft2 (for zero padding)

   % Normalization for correct boundary conditions
   % Center the (L1, L2) size mask within the maximal (L1max, L2max) extents.
   NP = real(fftshift(ifft2(abs(fft2(mask, L1max, L2max)).^2)));
   %NP = real(fftshift(ifft2(abs(fft2(mask, L1, L2)).^2)));
   % Collect the image areas (C1_A), image intensities (C1_N1 and C1_N2), FFT
   % transforms (C1_FF) and normalizations (C1_NP) separately, the last two
   % centered within the maximal rectangular extents (L1max, L2max), computing
   % the averaged C1 after the end of the loop.
   if corr_type == 'C'
      C1_A  = A;
      C1_N1 = N1;
      C1_N2 = N2;
      C1_FF = real(fftshift(ifft2(fft2(I1, L1max, L2max).* ...
                             conj(fft2(I2, L1max, L2max)))));
      C1_NP = NP;
      %C1 = A^2/N1/N2*real(fftshift(ifft2(fft2(I1, L1, L2).* ...
      %                              conj(fft2(I2, L1, L2)))))./NP;
   else
      % 2D G(r) with proper normalization
      C1_A  = A;
      C1_N1 = N1;
      C1_N2 = N1;
      C1_FF = real(fftshift(ifft2(abs(fft2(I1, L1max, L2max)).^2)));
      C1_NP = NP;
      %C1 = A^2/N1^2*real(fftshift(ifft2(abs(fft2(I1, L1, L2)).^2)))./NP;
   end

   sum_A  = sum_A  + C1_A ;
   sum_N1 = sum_N1 + C1_N1;
   sum_N2 = sum_N2 + C1_N2;
   sum_FF = sum_FF + C1_FF;
   sum_NP = sum_NP + C1_NP;
end
C1 = sum_A^2 / (sum_N1 * sum_N2) * sum_FF ./ sum_NP;

%only return valid part of G:
% a square with sides 2*rmax+1 about the center pixel, that is, the minimal
% rectangular region corresponding to the overlap of all the ROIs.
C = imcrop(C1, [floor(L2max/2+1)-rmax, floor(L1max/2+1)-rmax, 2*rmax, 2*rmax]);
%C = imcrop(C1, [floor(L2/2+1)-rmax, floor(L1/2+1)-rmax, 2*rmax, 2*rmax]);

xvals = ones(1, 2*rmax+1)'*(-rmax:rmax); %map to x positions with center x=0
yvals = (-rmax:rmax)'*ones(1, 2*rmax+1); %map to y positions with center y=0
zvals = C;

% convert x, y to polar coordinates
[theta,r,v] = cart2pol(xvals, yvals, zvals);

% Label each pixel in the minimal square region by its radius r from the
% center, then sort the r's and bin them.  Collect the pixels with radii in
% each bin and average the results to compute c(r).
Ar = reshape(r,1, (2*rmax+1)^2);
Avals = reshape(v,1, (2*rmax+1)^2);
[rr,ind] = sort(Ar);                  % sort by r values
vv = Avals(ind);                      % reindex g
r = 0:floor(max(rr));                 % the radii you want to extract
[n bin] = histc(rr, r-.5);            % bin by radius

for j = 1:rmax+1;                     % now get averages
   m = bin==j;
   n2 = sum(m);                       % the number of pixels in that bin
   if n2==0, vals(j)=0; er(j)=0;      % if no bins, no data
   else
      c(j) = sum(m.*vv)/n2;           % the average G values in this bin
                                      % the variance of the mean
      dc(j) = sqrt(sum(m.*(vv-c(j)).^2))/n2;
   end
end

r = 0:rmax;

end
