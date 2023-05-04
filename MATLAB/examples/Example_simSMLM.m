%% Example of generating sythetic SMLM data

SaveDir = fullfile(tempdir, 'smite', 'examples', 'simSMLM');
if ~isfolder(SaveDir)
   mkdir(fullfile(tempdir, 'smite'));
   mkdir(fullfile(tempdir, 'smite', 'examples'));
   mkdir(fullfile(tempdir, 'smite', 'examples', 'simSMLM'));
end

%Create sim object

S=smi_sim.SimSMLM()
NWings=20
S.NDatasets=20
S.SZ = 64;
S.simStar(NWings)

% Generate Images 
[Model,Data]=S.genImageStack();

% Generate Noisy Coordinates
[SMD_Noisy]=S.genNoisySMD()
figure;scatter(SMD_Noisy.X,SMD_Noisy.Y)

saveas(gcf, fullfile(SaveDir, 'sS1.png'));
