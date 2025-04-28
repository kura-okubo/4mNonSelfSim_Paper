# Paper for dynamics of non-self-similar earthquakes

This repository contains the notebooks for processing the experiments dataset, input files, and intermediate outputs for the paper of **"Dynamics of non-self-similar earthquakes illuminated by a controlled fault asperity"** (submitted).

[**Recipe of plotting figures in the article**](./plot_figures_recipe.md) is documented to reproduce the figures in the manuscript.


To reduce the repository size, we maintain the software used in this study across separate GitHub repositories as follows:

- [4mNonSelfSim_OpenSWPC](https://github.com/kura-okubo/4mNonSelfSim_OpenSWPC): modeling the waveform propagation through the four-meter-length rock specimens, by extending the functions of [OpenSWPC](https://openswpc.github.io) ([Maeda et al., 2017](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-017-0687-2)).

- [4mNonSelfSim_UGUCA](https://github.com/kura-okubo/4mNonSelfSim_UGUCA): performing dynamic rupture simulation of the non-self-similar laboratory earthquakes, adding the interface law to the SBIEM-based software [UGUCA](https://uguca.gitlab.io/uguca/) ([Kammer et al., 2019](https://www.sciencedirect.com/science/article/pii/S2352711021000959)).


## Dataset
The original dataset of stick-slip experiments is uploaded in Zenodo (https://doi.org/10.5281/zenodo.15233278), which includes high-sampling AE waveforms, strain, slip, and macroscopic measurement. We saved the data in MATLAB v7.3 format, extracting a 200 ms time window for each stick-slip main event from the continuous recording.

### Contents 
| Variable name | Description | Unit | Channels |
| --- | --- | --- | --- |
| `AEdatmat` | AE waveforms  | V (range: -10V to 10 V) | 32 |
| `AEsensor_x` | Sensor coordinate of AE sensor along fault | mm | | 
| `DXeast` | Macroscopic displacement measurement at east of lower rock specimen  | mm | 1 |
| `DXwest` | Macroscopic displacement measurement at west of lower rock specimen  | mm | 1 |
| `Disp_x` | Sensor coordinate of gap sensor along fault | mm | |
| `Dmat_event` | Slip measured by gap sensor | mm | 16 | 
| `Fs_AE` | Sampling frequency of AE sensor (10 MHz) | Hz | |
| `Fs_slip` | Sampling frequency of gap sensor (50 kHz) | Hz | |
| `Fs_strain` | Sampling frequency of strain gouge (1 MHz) | Hz | |
| `NPmacro` | Normal pressure measurements at flat jacks^1 | MPa | 8 |
| `SGB_x` | Sensor coordinate of biaxial strain gouge along fault | mm | |
| `SGT_x` | Sensor coordinate of biaxial strain gouge along fault | mm | |
| `SSmacro` | Macroscopic shear stress measured at shear jack | MPa | | 
| `Snmat` | Normal stress measured by strain gouge^2 | MPa | 32 |
| `Spmat` | Horizontal stress ($\\sigma_{xx}$) measured by strain gouge^2 | MPa | 32 |
| `Taumat2` | Shear stress change measured by biaxial strain gouge^2 | MPa | 32 |
| `Taumat3` | Shear stress change measured by triaxial strain gouge^2 | MPa | 32 |
| `Tstart` | Absolute start time of this event data  | sec | |
| `event_winlen` | Event window length of this data  | sec | |
| `tmat_AE_event` | Time vector of AE | sec | |
| `tmat_macro` | Time vector of macroscopic data | sec | |
| `tmat_slip_event` | Time vector of slip  | sec | |
| `tmat_strain_event` | Time vector of strain  | sec | |

- ^1 The normal pressure applied by the flat jacks is converted to normal stress on the fault by a factor of 2/3, accounting for the contact area between the flat jacks and the rock specimen.

- ^2 Strain measurements are calibrated at the start of the experiment when the top rock block is lifted, then the recording is initiated.
 
### How to load

In matlab, load the data for the stick-slip, for example, event id 29 by `A=load("eventdata_FB03_087_event29.mat")`.

## Reference
Okubo, K., Yamashita, F., & Fukuyama, E. (2025) Dynamics of non-self-similar earthquakes illuminated by a controlled fault asperity, submitted.
