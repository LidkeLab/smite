%Example4
% -- needs the matlab toolbox dipimage, free download at www.diplib.org
% run each section separtely. This file shows how the functions "imres_ims"
% and "imres_locs" can be used to generate and automatically display results
% from positions or images.
% -- you need to have the fullpaht to the FIREfunctions directory on the MATLAB path  

% this examples takes a few minutes to complete

clear all
close all

%% loading the data
% localizations of single emitters in the format:
%  x, y, t [pixels, pixels, frames]

coords = dlmread(['..' filesep 'ExampleData' filesep 'example_Fig2a.dat'],',');
pixelsize = 16e3/150;                                                           %in nanometers

%% Select input format and outputs
positions_order = 'xyt';
positions_units = 'CCD pixels';
show_im = 1;
show_frc = 1;
show_timefractions = 1;

%% call imres_locs
superzoom = 10;
sz = superzoom * 256 * [1 1];
SR_pixelsize = pixelsize/superzoom;
nblocks = 20;
timefractions = 10;
reps = 20;

[res_value frc_curve im_out res_high res_low res_t] = imres_locs(coords,sz,superzoom,nblocks,timefractions,reps,SR_pixelsize,positions_order, positions_units, show_im, show_frc,show_timefractions);
% [res_value frc_curve im_out res_high res_low res_t] = imres_locs(coords,sz,superzoom,nblocks,timefractions,reps,SR_pixelsize);

%% call imres_ims
maxt = max(coords(:,3));

% Define two images for resolution computation
in1 = binlocalizations(coords(coords(:,3)<maxt/2,:),sz(1),sz(2),superzoom);
in2 = binlocalizations(coords(coords(:,3)>=maxt/2,:),sz(1),sz(2),superzoom);

% Use res_ims
[res_value frc_curve res_high res_low] = imres_ims(in1,in2,SR_pixelsize,show_frc);
