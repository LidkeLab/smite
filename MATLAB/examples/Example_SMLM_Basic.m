%Demonstrate basic functionality of SMLM Analysis

%Create a SMLM object
SMLMobj=smi.SMLM()  %no input pops open a gui

%Use gui to navigate to a test dataset such as this TIRF DNA-PAINT: 

%  Y:\Sandeep\20-11-2020-DNA_PAINT_Tubulin\Dock2-Cell1-2020-11-12-10-29-58.h5

% Set other SMF parameters and try a test fit from gui



%% Scripted Analysis

%or set it manually:
SMLMobj.SMF.Data.FileDir='Y:\Sandeep\20-11-2020-DNA_PAINT_Tubulin';
SMLMobj.SMF.Data.FileName='Dock2-Cell1-2020-11-12-10-29-58.h5';

%give a Analysis ID
SMLMobj.SMF.Data.FileName='Dock2-Cell1-2020-11-12-10-29-58.h5';




