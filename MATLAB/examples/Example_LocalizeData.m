%% Demonstrates LocalizeData Class

%create a test dataset
B = smi_sim.GaussBlobs.genRandomBlobImage();
B = poissrnd(B);

%setup SMF
SMF = smi_core.SingleMoleculeFitting()

%setup LocalizeData object
LD = smi_core.LocalizeData(B, SMF)

%Localize with defaults
[SMD] = LD.genLocalizations();

%Set Verbose to give color overlay output
LD.Verbose = 3;
[SMD] = LD.genLocalizations();
