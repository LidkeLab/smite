%Example3: anisotropy extention in 3D
% Fourier Plane Correlation
% -- needs the matlab toolbox dipimage, free download at www.diplib.org

clear all
close all

%% loading the data
% localizations of single emitters in the format:
%  x, y, t [pixels, pixels, frame]
coords = dlmread(['..' filesep 'ExampleData' filesep 'example_Fig4.dat'],',');
pixelsize = 16e3/147;                                                           %in nanometers xy, in z microns
SRpixelsize = 10;


%% make two images
zoom = 1;
coordsin = coords./SRpixelsize;
numpixelsx = ceil(max(coordsin(:,1)));
numpixelsy = ceil(max(coordsin(:,2)));
numpixelsz = ceil(max(coordsin(:,3)));
Nem = length(coordsin);
im1 = binlocalizations3D(coordsin(1:round(Nem/2),1:3),numpixelsx,numpixelsy,numpixelsz,zoom);
im1 = mat2im(im1);
im2 = binlocalizations3D(coordsin(round(Nem/2)+1:end,1:3),numpixelsx,numpixelsy,numpixelsz,zoom);
im2 = mat2im(im2);
Nim1 = imsize(im1);

% crop images
Nsize = 512;
xmin = 90;
ymin = 0;
zmin = 0;
imcrop1 = cut(im1,Nsize,[xmin ymin zmin]);
imcrop2 = cut(im2,Nsize,[xmin ymin zmin]);
szCR = imsize(imcrop1);
imcrop1 = extend(imcrop1,[szCR(1) szCR(1) szCR(1)]);
imcrop2 = extend(imcrop2,[szCR(1) szCR(1) szCR(1)]);
clear im1 im2

%% compute FPC slices only

[fpc_xy,fpc_xz,fpc_yz]=fpc(imcrop1,imcrop2);
dipshow(fpc_xy,'lin')
dipshow(fpc_xz,'lin')
dipshow(fpc_yz,'lin')


%% compute the full FPC (much slowers as the full 3D volume is computed)
% this routine requires much memory if the images are large, you need to
% Fourier transform the cubes NxNxN and operate on them
% for the manuscript a much larger part of the image was used and more
% angles, this is very time costly.

[~, ~, ~, fpc_out] = fpc(imcrop1, imcrop2, 360*3, 0.2);
sz=imsize(fpc_out);
h=dipshow(fpc_out,'lin');
dipmapping(h,'global','slice',round(sz(1)/2))

