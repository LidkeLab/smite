load('C:\Users\David\Documents\GitHub\smite\MATLAB\+smi_stat\testTR.mat')

%%
DE = smi_stat.DiffusionEstimator(TR);
DE.FitMethod = 'WeightedLS';
DE.UnitFlag = 0;
% DE.FitTarget = 'CDFOfJumps';
DE.FitTarget = 'MSD';
DE.estimateDiffusionConstant();
DE.plotEnsembleMSD([], DE.MSDEnsemble, DE.DiffusionStruct, DE.UnitFlag)



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