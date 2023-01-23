PC = smi_cluster.PairCorrelation();
PC.ResultsDir = 'Results';

% Generated data.
SZ = 1000;
x1 = SZ * rand(2 * SZ, 1);
y1 = SZ * rand(2 * SZ, 1);
XY1 = [x1, y1];
SMD1.X = x1;
SMD1.Y = y1;
x2 = [SZ * rand(3 * SZ, 1)];
y2 = [SZ * rand(3 * SZ, 1)];
XY2 = [x2, y2];
SMD2.X = x2;
SMD2.Y = y2;

PC.ROI = [0, SZ, 0, SZ];   % [x_min, x_max, y_min, y_max]
% These two numbers below are often the same, but the user can increase
% HistBinSize to make bigger internal pixels (or bigger internal histogram
% image bins) and so produce more smoothing, or decrease this quantity and
% attempt greater detail.
PC.PixelSize = 1;         % Actual camera pixel size (nm).
PC.HistBinSize = 1;       % Internal image pixel size (nm).

% Make a RoI structure (see also smi_cluster.ROITools).  This is typically
% used for invoking pair_correlation_ROIcombined with multiple ROIs, which
% are combined.  The size of each ROI need not be the same, although
% cleaner results will be produced if they are all the same size.
n_ROIs = 1;
RoI{1}.ROI = PC.ROI;
RoI{1}.X   = {XY1(:, 1), XY2(:, 1)};
RoI{1}.Y   = {XY1(:, 2), XY2(:, 2)};

% Typically, pair_correlation is used for comparing two images/ROIs, while
% pair_correlation_ROIcombined is used when combining the pair correlation
% results for several pairs of images/ROIs.  pair_correlation_Veatch is
% basically the original Sarah Veatch code left for comparison purposes.
% The results of all three routines in these examples should produce
% similar results (identical for the first two when n_ROIs = 1).
% NOTE: SMD structures or XY matrices can be provided as input.
PC.BaseName = '9021';
results_pcc  = PC.pair_correlation(SMD1, SMD2)
results_pcc  = PC.pair_correlation(XY1, XY2)
results_Rpcc = PC.pair_correlation_ROIcombined(2, n_ROIs, RoI)
results_Vpcc = PC.pair_correlation_Veatch(SMD1, SMD2, 'cross')
results_Vpcc = PC.pair_correlation_Veatch(XY1, XY2, 'cross')

PC.BaseName = '9021_5';
results_pac1  = PC.pair_correlation(SMD1)
results_pac1  = PC.pair_correlation(XY1)
results_Rpacc = PC.pair_correlation_ROIcombined(1, n_ROIs, RoI, 1)
results_Vpac1 = PC.pair_correlation_Veatch(SMD1,  [], 'auto')
results_Vpac1 = PC.pair_correlation_Veatch(XY1,  [], 'auto')

PC.BaseName = '9021_10';
results_pac2  = PC.pair_correlation(SMD2)
results_pac2  = PC.pair_correlation(XY2)
results_Rpac2 = PC.pair_correlation_ROIcombined(1, n_ROIs, RoI, 2)
results_Vpac2 = PC.pair_correlation_Veatch(SMD2, [], 'auto')
results_Vpac2 = PC.pair_correlation_Veatch(XY2, [], 'auto')

fprintf('Done pair correlation.\n');
fprintf('All results in %s\n', PC.ResultsDir);
