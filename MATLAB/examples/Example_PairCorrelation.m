PC = smi_cluster.PairCorrelation();
PC.ResultsDir = 'Results';

% Data from some external source.
[x1, y1] = textread('9021_5.txt',  '%*u %u %u %*u', 'headerlines', 1);
XY1 = [x1, y1];
[x2, y2] = textread('9021_10.txt', '%*u %u %u %*u', 'headerlines', 1);
XY2 = [x2, y2];

PC.ROI = [0, 7400, 0, 6000];   % [x_min, x_max, y_min, y_max]
% These two numbers below are often the same, but the user can increase
% HistBinSize to make bigger internal pixels (or bigger internal histogram
% image bins) and so produce more smoothing, or decrease this quantity and
% attempt greater detail.
PC.PixelSize = 2.7559;         % Actual camera pixel size (nm).
PC.HistBinSize = 2.7559;       % Internal image pixel size (nm).

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
PC.BaseName = '9021';
results_pcc  = PC.pair_correlation(XY1, XY2)
results_Rpcc = PC.pair_correlation_ROIcombined(2, n_ROIs, RoI)
results_Vpcc = PC.pair_correlation_Veatch(XY1, XY2, 'cross')

PC.BaseName = '9021_5';
results_pac1  = PC.pair_correlation(XY1)
results_Rpacc = PC.pair_correlation_ROIcombined(1, n_ROIs, RoI, 1)
results_Vpac1 = PC.pair_correlation_Veatch(XY1,  [], 'auto')

PC.BaseName = '9021_10';
results_pac2  = PC.pair_correlation(XY2)
results_Rpac2 = PC.pair_correlation_ROIcombined(1, n_ROIs, RoI, 2)
results_Vpac2 = PC.pair_correlation_Veatch(XY2, [], 'auto')