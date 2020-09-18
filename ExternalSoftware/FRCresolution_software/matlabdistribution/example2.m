%Example2: anisotropy extention in 2D
% Fourier Line Correlation
% -- needs the matlab toolbox dipimage, free download at www.diplib.org

clear all
close all

%% loading the data
% localizations of single emitters in the format:
%  x, y, t [pixels, pixels, frame]
coords = dlmread(['..' filesep 'ExampleData' filesep 'example_Fig2j.dat'],',');     %1/10 of the total data
pixelsize = 16e3/160;                                                               %in nanometers

superzoom = 5;                                                                      %10 was used in the paper
szx = superzoom * 180;
szy = superzoom * 180;

%%make two images
t = coords(:,3);
tmax = max(t);
in1 = binlocalizations(coords(t<tmax/2,:), szx, szy, superzoom);
in2 = binlocalizations(coords(t>=tmax/2,:), szx, szy, superzoom);

%%compute the FLC, (takes ~30 seconds)
flc_out = flc(in1,in2);

flc_cut = cut(flc_out, round(imsize(flc_out)./4));
h = dipshow(flc_cut,'lin');
dipmapping(h,'colormap',jet)


%% show via matlab toolbox
allfreqs = linspace(-1/8,1/8,imsize(flc_cut,1))/(pixelsize/superzoom);
[qx, qy] = meshgrid(allfreqs,allfreqs);
figure;
pcolor(qx,qy, im2mat(rot90(flc_cut)))
xlabel('q_x (nm^{-1})')
ylabel('q_y (nm^{-1})')
shading flat;
colormap('jet');
axis equal
%colorbar;
