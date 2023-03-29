### CalibrationFile

Example camera calibration file:

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
