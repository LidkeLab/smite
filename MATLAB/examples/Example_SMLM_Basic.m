%Demonstrate basic functionality of SMLM Analysis

%Create a SMLM object
SMLMobj=smi.SMLM()  %no input pops open a gui

%Use gui to navigate to a test dataset, for example a TIRF DNA-PAINT .h5 file.

% Set other SMF parameters and try a test fit from gui



%% Scripted Analysis

%or set it manually:
%SMLMobj.SMF.Data.FileDir='/path/to/data';
%SMLMobj.SMF.Data.FileName='Data.h5';

%give a Analysis ID
%SMLMobj.SMF.Data.AnalysisID='ID';
