# Recipe of plotting the figures

!!! note "Notebooks and codes to reproduce the figures"
    This page shows the list of the link to the notebooks and codes to plot the figures. We cleaned up the output of the cells showing the figures in the notebooks due to minimizing the size of repositories and the copyrights. To replot them, please download the data and follow the instructions documented in the directories.


## Figure Reproduction Notes

- **Figure 1a.** Conceptual illustration of non-self-similar earthquakes.  
  Run [`01_nonselfsimilar_schematic.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/Others/Fig1_introduction/Fig1a_nonselfsimilar_schematic/01_nonselfsimilar_schematic.ipynb)

- **Figures 1b–c.** Diagram of the four-meter biaxial friction apparatus.  
  Run [`01_Figs1bc_4mbiaxapparatus_schematic.ipynb`](/Users/kokubo/Dropbox/NIED_RESEARCH/4mBIAX_submission/4mNonSelfSim_Paper_v2/Others/Fig1_introduction/Fig1b_4mbiaxapparatus_schematic/code/01_Figs1bc_4mbiaxapparatus_schematic.ipynb)

- **Figures 1d–f.** AE waveforms of GP events with zoomed-in windows.  
  Run [`01_Introduction_plot_gougeevents.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/ComputeScaling/code/01_Introduction_plot_gougeevents.ipynb)

- **Figure 2a.** Shear stress and AE waveforms for a representative GP event.  
  Run [`plot_example_master_fb03_087_event29_Fig2_master.m`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/PlotEvent/code/plot_example_master_fb03_087_event29_Fig2_master.m)

- **Figure 2b.** Cumulative slip during the GP event.  
  Run [`plot_eventcumulativeslip_event29_forFig2.m`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/PlotEvent/code/plot_eventcumulativeslip_event29_forFig2.m)

- **Figure 2c.** Full waveform view of a GP event.  
  Run [`Fig2_plot_eventwaveform.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/Others/Fig2_plot_eventwaveform/code/Fig2_plot_eventwaveform.ipynb)

- **Figures 3a–b.** Source time function (STF) fitting for GP events.  
  Run [`06_plotfittingSTF.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/ComputeScaling/code/06_plotfittingSTF.ipynb)

- **Figure 3c.** Moment-duration scaling plot.  
  Run [`07_loglinearfit_STF.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/ComputeScaling/code/07_loglinearfit_STF.ipynb)

- **Figures 4a-b.** Initial stress distribution in the rupture simulation.  
  Run [`02_plot_initialcondition.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/RuptureSimulation/main_casestudy/postprocess_dynrup/code/02_plot_initialcondition.ipynb)

- **Figure 4c.** Shear traction history during rupture.  
  Run [`04_plot_sheartractionhistory.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/RuptureSimulation/main_casestudy/postprocess_dynrup/code/04_plot_sheartractionhistory.ipynb)

- **Figures 4d–f.** Snapshots of rupture front propagation.  
  Run [`03_plot_snapshots_mastercase.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/RuptureSimulation/main_casestudy/postprocess_dynrup/code/03_plot_snapshots_mastercase.ipynb)

- **Figure 5.** STF of dynamic rupture simulations.  
  Run [`05_plot_masterSTF.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/RuptureSimulation/main_casestudy/postprocess_dynrup/code/05_plot_masterSTF.ipynb)

---

## Supplementary Figure Reproduction Notes

- **Figure S1.** Schematic of the experimental frame and AE sensor array.  
  Run [`plot_4mframeandSensors.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/Others/SensorArray/code/plot_4mframeandSensors.ipynb)

- **Figure S2.** Photographs of the experimental apparatus and gouge patch (GP).  
  *No script required.*

- **Figure S3.** Histogram of GP event counts and spatial distribution of normal stress.  
  Run [`plot_GPeventActivity_and_NormalStress.m`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/Others/GPeventsActivity/code/plot_GPeventActivity_and_NormalStress.m)

- **Figure S4.** AE waveforms of GP events across all AE sensors.  
  Run [`01_Introduction_plot_gougeevents.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/ComputeScaling/code/01_Introduction_plot_gougeevents.ipynb)

- **Figure S5.** Shear stress and AE waveform of a GP aftershock event.  
  Run [`plot_example_master_fb03_087_event35_FigSupp.m`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/PlotEvent/code/plot_example_master_fb03_087_event35_FigSupp.m)

- **Figure S6.** Shear stress and AE waveform of another GP aftershock event.  
  Run [`plot_example_master_fb03_087_event35_FigSupp.m`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/PlotEvent/code/plot_example_master_fb03_087_event35_FigSupp.m)

- **Figures S6a–b.** Schematic diagrams for estimating attenuation factors.  
  Run [`08_compute_attenuation_factor.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/Calibration/Attenuation/code/08_compute_attenuation_factor.ipynb)

- **Figure S6c.** Quantitative estimation of attenuation (1/Q).  
  Run [`09_compute_Qinv.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/Calibration/Attenuation/code/09_compute_Qinv.ipynb)

- **Figure S7.** Spectral ratio analysis of GP events.  
  Run [`02_stacked_spectralratioanalysis_nonselfsimilar.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/Others/SpectralRatio/code/02_stacked_spectralratioanalysis_nonselfsimilar.ipynb)

- **Figure S8.** Cumulative local slip distributions across the fault.  
  Run [`05_localslip_slipvel_masterplot.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/GougeEventStats/M0_LocalSlip_and_SlipVel/code/05_localslip_slipvel_masterplot.ipynb)

- **Figure S9.** STF of dynamic rupture simulation without self-healing.  
  Run [`05_plot_masterSTF.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/RuptureSimulation/main_casestudy/postprocess_dynrup/code/05_plot_masterSTF.ipynb)

- **Figure S10.** Statistical summary of gouge patch activity.  
  See [`GougePatch`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/Experiments/GougePatch)

- **Figure S11.** Schematic of AE sensor calibration setup.  
  *No script required.*

- **Figure S12.** Flowchart of AE sensor calibration procedure.  
  *No script required.*

- **Figure S13.** AIC-based model selection in AE sensor calibration.  
  Run [`plot_modelAIC.m`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/Calibration/AEsensor_Calibration/Evaluation_PAZmodel/code/plot_modelAIC.m)

- **Figure S14a.** Bode plot of transfer function.  
  See [`AEsensor_Calibration_ARX`](https://github.com/kura-okubo/AEsensor_Calibration_ARX)

- **Figure S14b.** Comparison of responses before and after correction.  
  Run [`validate_metalblock_PAZmodel.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/Calibration/AEsensor_Calibration/Evaluation_PAZmodel/code/validate_metalblock_PAZmodel.ipynb)

- **Figure S15.** Locations of ball-drop calibration tests.  
  Run [`02_plot_and_makelocationtable_v01.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/Calibration/SensorCoupling_BallDrop/code/02_plot_and_makelocationtable_v01.ipynb)

- **Figure S16.** Schematic of waveform modeling and calibration strategy.  
  See [`4mNonSelfSim_OpenSWPC – cross-verification`](https://github.com/kura-okubo/4mNonSelfSim_OpenSWPC/tree/develop/cross-verification)

- **Figure S17.** Results of calibration for waveform propagation modeling.  
  See [`4mNonSelfSim_OpenSWPC – cross-verification`](https://github.com/kura-okubo/4mNonSelfSim_OpenSWPC/tree/develop/cross-verification)

- **Figure S18.** Waveform fitting for sensor coupling calibration.  
  Run [`08_plot_surfaceeffect_result_Case2.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/Calibration/SensorCoupling_BallDrop/code/08_plot_surfaceeffect_result_Case2.ipynb)

- **Figure S19.** Schematic illustration of the aperture effect.  
  *No script required.*

- **Figure S20.** Sensor coupling calibration results.  
  Run [`08_plot_surfaceeffect_result_Case2.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/Calibration/SensorCoupling_BallDrop/code/08_plot_surfaceeffect_result_Case2.ipynb)

- **Figure S21.** STF fitting results for all GP events.  
  Run [`06_plotfittingSTF.ipynb`](https://github.com/kura-okubo/4mNonSelfSim_Paper/blob/dev/ComputeScaling/code/06_plotfittingSTF.ipynb)

---