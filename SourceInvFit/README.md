# Documentation of Moment Tensor inversion tools

## AE detection and location
We made a set of scripts and notebooks to evaluate the seismic moment tensor and the duration of gouge-mediated seismic events.

### 01_AE_locateevent_GUI.py
Pick the first arrival time to determine the event locations using the GUI.
We picked around four stations, and conducted the grid search to find the most likely location. We checked if it consistents to the location of gouge patches, and saved the picked time into `data/AElocation/arrivalpick`.

### 02_AE_relocation.ipynb
To better fit between the data and numerical waveforms, we visually adjusted the first arrival time by relocating the source locations. The output is stored in `data/AElocation/relocation`. `02_AE_relocation_all.py` is the script to run for all the gouge events.

### 03_AE_save_eventwaveform_relocated.ipynb
We read the AE raw waveform from the event `.mat` data with the origin time determined with the relocation. We stored the set of AE waveform into the obspy stream, and dumped the data in the `data/03_AEobs_waveform`.　We also dumped the location datasheet into `data/datacsv/AE_obs_location.csv`

### 04_AE_convert_to_isoparametric_coordinate.ipynb
We shifted the coordinates of source location such that we conduct the numerical simulation of wave propagation on the isoparametric coordinate system. The input file of OpenSWPC green function's table is output in `data/datacsv/green_in_AEevent_biax.txt`.

### 05_numericalsimulation.ipynb
We conducted the numerical simulation using OpenSWPC. We set the Green's function mode in the input file, and run in the HPC cluster. The Green's function is in `out/green/S00` in the sac format.

### 06_assemble_greensfunction_MTinv_removeresp.ipynb
We gather the observed waveform with the correction of instrumental response and the sensor coupling by the ball-drop test and the synthetic Green's function associated with the events. The stream is output into `"data/06_assemble_gf`. 

### 07_AEevents_GridSearchMTinversion_v09_PS_removeresp.ipynb
We conducted the source inversion by fitting the waveform. `07_2_allevents_AEevents_GridSearchMTinversion_v09_PS_removeresp.py` is the script version of the notebook to run the loop of all of the events. The result is output in `data/07_DATA_MTinversion`.

### 08_summarize_gougeevent_stats.ipynb
We plotted the contour of VR to check the convergence of source inversion. Then, output `data/datacsv/gridsearch_bestparam_M0andTR_fb03-087.csv`, which includes the best fit parameters associated with $M_0$ and $T_R$. Note that the seismic moment is NOT $\hat{M}_0$. The magnitude is scated as $M_0 = \hat{M}_0/\sqrt{2}$.







