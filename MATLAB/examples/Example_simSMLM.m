
%% Example of generating sythetic SMLM data

%Create sim object

S=smi_sim.SimSMLM()
NWings=20
S.NDatasets=20
S.simStar(NWings)

% Generate Images 
[Model,Data]=S.genImageStack();

% Generate Noisy Coordinates
[SMD_Noisy]=S.genNoisySMD()
figure;scatter(SMD_Noisy.X,SMD_Noisy.Y)




