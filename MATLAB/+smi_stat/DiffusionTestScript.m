load('C:\Users\David\Documents\GitHub\smite\MATLAB\+smi_stat\testTR.mat')

%%
DE = smi_stat.DiffusionEstimator(TR);
DE.FitMethod = 'WeightedLS';
DE.UnitFlag = 1;
DE.MaxFrameLag = 5;
DE.FitTarget = 'CDFOfJumps';
% DE.FitTarget = 'MSD';
DE.FitIndividualTrajectories = false;
DE.Verbose = 1;
DE.estimateDiffusionConstant();

figure();
DE.plotEnsembleMSD([], DE.MSDEnsemble, DE.DiffusionStruct, ...
    DE.DiffusionModel, DE.UnitFlag)

DE.plotEnsembleCDFOfJumps([], DE.MSDEnsemble, DE.DiffusionStruct, ...
    DE.DiffusionModel, DE.UnitFlag)



[MSDSingleTraj, MSDEnsemble] = ...
    DE.computeMSD(DE.TR, DE.MaxFrameLag, DE.Verbose);
Weights = MSDEnsemble.NPoints ./ MSDEnsemble.FrameLags.^2;
[BetaHat, BetaHatSE] = smi_stat.leastSquaresFit(...
    MSDEnsemble.FrameLags, MSDEnsemble.MSD, Weights)

FitResults = fit(MSDEnsemble.FrameLags, MSDEnsemble.MSD, 'poly1', ...
    'Weights', Weights);
FitParams = coeffvalues(FitResults);
FitParamsCI = confint(FitResults, erf(1 / sqrt(2)));
FitParamsSE = FitParamsCI(2, :) - FitParams

%%
X = linspace(0, 100, 1e3).';
Y = 1.23 + 3.21*X + (0.1^2)*randn(1e3, 1);
[ParamsHat, ParamsHatSE] = smi_stat.bootstrapFit(X, Y, [1, 1]);

