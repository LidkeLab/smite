### CalibrationFile

For quantitative analysis of single molecule localization microscopy
(SMLM) images, pixelwise properties of the camera must be well
characterized and accounted for in the analysis algorithm to avoid
artifacts. Offset and gain maps are considered the basis for numerous
camera correction algorithms. We collect the offset and gain maps
of sCMOS cameras and add them via the ***smite*** GUI or in an analysis
script by providing the path of a saved `CalibrationFile`.
We determine the mean ('offset'),
variance and amplification gain of each pixel of the signal over many
images at several light levels. Following are the necessary steps where
we use
[***MIC***](https://github.com/LidkeLab/matlab-instrument-control.git)),
[DIPimage](https://diplib.org/) and ***smite*** tools:

* Connect to the devices:

    1. Connect to camera:
```
    CameraSCMOS = MIC_HamamatsuCamera();
    CameraSCMOS.ReturnType = 'matlab';
    CameraSCMOS.gui();
```
    2. Connect to lamp:
```
    Lamp660 = MIC_ThorlabsLED('Dev1', 'ao0');
    Lamp660.gui();
```
* Define a few parameters: 
```
    LampPowerRange = linspace(0, 3.6, 20); % selected to reach ~600 max. camera counts
    CameraSCMOS.DefectCorrection = 1; % no correction
    CameraSCMOS.ScanMode = 1; % slow scan
    CameraSCMOS.ExpTime_Sequence = 0.01;
    CameraSCMOS.SequenceLength = 1000;
    CameraSCMOS.ROI = [897, 1152, 897, 1152]; % [XStart, XEnd, YStart, YEnd]
```
* Collect the gain/offset data:
```
    Params = [];
    MeanLevel = [];
    VarLevel = [];
    CameraSCMOS.AcquisitionType = 'sequence';
    CameraSCMOS.ExpTime_Sequence = 0.01;
    CameraSCMOS.setup_acquisition();
    for ii = 1:numel(LampPowerRange)
       fprintf('Lamp power %i of %i\n', ii, numel   (LampPowerRange))
       Lamp660.setPower(LampPowerRange(ii));
       pause(1)
       CameraSCMOS.start_sequence();
       MeanLevel = cat(3, MeanLevel, mean(CameraSCMOS.Data, 3));
       VarLevel = cat(3, VarLevel, var(single(CameraSCMOS.Data), [], 3));
    end
    Lamp660.setPower(0);
    Params.MeanLevel = single(MeanLevel);
    Params.VarLevel = single(VarLevel);
```
* Perform Least-Squares Fit on collected data:

    leastSquaresFit from smi_stat class performs a least squares (ls)
    fit on the provided data.
    This function computes a least squares fit for the paired data in
    XData, YData, and Weights. Check
    [`smi_stat.leastSquaresFit()`](../../MATLAB/+smi_stat/leastSquaresFit.m)
    for further documentation.
```
    % Here [:, :, 1] is ls offset, [:, :, 2] is ls slope (Gain)
    Beta = NaN(size(MeanLevel, 1), size(MeanLevel, 2), 2);
    for ii = 1:size(MeanLevel, 1)
       disp(ii)
       for jj = 1:size(MeanLevel, 2)
          Beta(ii, jj, 1:2) = smi_stat.leastSquaresFit(                   ...
             squeeze(MeanLevel(ii, jj, :)), squeeze(VarLevel(ii, jj, :)), ...
             1 ./ squeeze(VarLevel(ii, jj, :)));
       end
    end
    dipshow(Beta(:, :, 2))

    figure;
    histogram(Beta(:, :, 2))
    xlabel('Gain (ADU/e-)')

    figure;
    histogram(MeanLevel(:,:,1))
    xlabel('ADU')
    Params.CCDVar=single(VarLevel(:,:,1));
    Params.Gain=single(Beta(:,:,2));
    Params.CCDOffset=single(MeanLevel(:,:,1));
```
* Save the calibration:
```
    SaveDir = 'Y:\sCMOS Calibrations\Sequential SR';
    FileName = fullfile(SaveDir, ...
       ['GainCalibration-', smi_helpers.genTimeString()]);
    Params.CameraObj.ROI = CameraSCMOS.ROI;
    Params.LampPowerRange = single(LampPowerRange);
    save(FileName, 'Params', '-v7.3')
```
* Example camera calibration file:

```
   >> CalibrationFile.Params

   ans = 

     struct with fields:

                      MeanLevel: [256×256×20 single]
                       VarLevel: [256×256×20 single]
                 DarkImagesMean: [256×256×5 single]
                  DarkImagesVar: [256×256×5 single]
              DarkImagesExpTime: [0.0100 0.1325 0.2550 0.3775 0.5000]
       DarkImagesSequenceLength: 1000
                         CCDVar: [256×256 single]
                        CCDGain: [256×256 single]
                      CCDOffset: [256×256 single]
                      CameraObj: [1×1 struct]
                 LampPowerRange: [1×20 single]
                           Gain: [256×256 single]
```
