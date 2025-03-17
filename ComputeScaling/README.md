# Compute moment-duration scaling from the P wave pulse of the gouge event waveform

## 01_Introduction_plot_gougeevents.ipynb
We repick the P wave onset using the second derivative of the waveform for the accurate analysis of the source parameters on the GP-generated events. We need to run this notebooks four times as `repeated_sensor=OL07`, `OL08`, `OL22`, `OL23`, to complete the process. 

## 02_trim_Pwave_window.ipynb
Trimming the P wave window after correcting for the aperture effect factor. We also need to run four times as `sensor_id=7, 8, 22, 23` to complete the process. 

## 03_compute_displacementpulse.ipynb
Computing the displacement pulse with preprocessing including the detrend, low-pass filtering, correction of attenuation, integration, and align the onset.The noise level is also evaluated, used for the threshold to evaluate the quality of source parameter estimation.

## 04_fittingSTF_all.ipynb
Fitting the synthetic cosine STF to the observed P wave displacement pulse.

## 05_compute_STFstats.ipynb
Computing the statistics of the source parameters.

## 06_plotfittingSTF.ipynb
Plot the result of fitting.

## 07_loglinearfit_STF.ipynb
Plot the moment-duration scaling of GP event. We performed the major axis regression to evaluate the scaling exponent.

