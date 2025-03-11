# Calibration of AE sensor with ARX model

## AEsensor_Calibration_DataConversion
We converted the raw binary record of the calibration waveforms associated with the LDV and AE sensors to the stacked waveforms, saved in the `Evaluation_PAZmodel/data`. We confirmed the reproducibility for those data on March 11, 2025.


## [AEsensor_Calibration_ARX](https://github.com/kura-okubo/AEsensor_Calibration_ARX)

We developed a script to evaluate the ARX model using AE sensor and LDV records, which is available in a dedicated repository [here](https://github.com/kura-okubo/AEsensor_Calibration_ARX).

## Evaluation_PAZmodel
This directory contains the optimization of poles and zeros based on the AIC, and the plots of results. We also checked the linearity of the piezoelectric source between 100V and 200V step pulse using LDV. The robustness was evaluated using three different measurement locations.
