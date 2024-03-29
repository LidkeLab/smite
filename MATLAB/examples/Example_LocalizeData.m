%% Demonstrates LocalizeData Class

SaveDir = smi_helpers.mkSMITETmpDir('examples', 'LocalizeData');

%create a test dataset
B = smi_sim.GaussBlobs.genRandomBlobImage();
B = poissrnd(B);

%setup SMF
SMF = smi_core.SingleMoleculeFitting()

%setup LocalizeData object
LD = smi_core.LocalizeData(B, SMF)

%Localize with defaults
[SMD] = LD.genLocalizations();

saveas(gcf, fullfile(SaveDir, 'LD1.png'));

%Set Verbose to give color overlay output
LD.Verbose = 3;
[SMD] = LD.genLocalizations();

saveas(gcf, fullfile(SaveDir, 'LD2.png'));
