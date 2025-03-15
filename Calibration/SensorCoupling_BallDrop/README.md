# Calibration of Sensor coupling using ball-drop impact

## Workflow

0. We conducted the ball drop test at 32 locations on the fault. `data/DATA_RAW` contains the raw data of the waveform associated with ball drop impact for 3 times at each location.

1. `01_relocate_balldrop_v01.ipynb` conducts the picking for the onset of the P wave and relocating the ball drop impact by the grid search. It dumps the stream including the locations and set of waveforms near the source. We removed the instrumental response at the beginning of processing in this stage. Therefore, the the gain can be obtained by the coupling factor and scale factor associated with the response removal. The data is output into `data/DATA_computeCF`.

2. `02_plot_and_makelocationtable_v01.ipynb` plots the results of ball drop test and generates the csv table for the relocated locations of ball impact (`balldrop_locations.csv`).

3. `03_convert_to_isoparametric_coordinate.ipynb` converts the coordinate to the isoparametric locations, i.e., relative coordinates from the source to AE sensor. It saves the traces to the pickle with the isoparametric coordinates to `DATA_isocoord`. The Green's function table as an input file for the numerical modeling with OpenSWPC is also output. `aux03_1_plot_isoloc_greensfunctionsource_all.ipynb` plots all the ball drop locations with the isoparametric coordinate.

4. `04_numericalmodeling_waveform` contains the input files to run the OpenSWPC with the Green's function mode. You can run it locally with `mpirun -np 16 swpc_3d.x -i ./input.inf`, or submit a job in the cluster. 

5. `05_assemble_greensfunction.ipynb` assembles the syntetic and observed waveforms with respect to the AE sensors. We preprocessed the numerically modeled waveforms; converted from displacement to the velocity by differenciation, and time shift with the amount of pretrigger with zero padding. `aux05_1_compare_gain_dispvel_to_previousstudies.ipynb` conducted the comparison on the sensitivity to the displacement to verify it with Wu and McLaskey (2018). `aux05_2_check_balldrop_sourcetimefunction.ipynb` plots the source time function obtained from the Hertzian impact theory.

6. `06_surfaceefect_assemble_Aij.ipynb` to compute the P wave amplitude for the synthetic and observed waveforms. `Case 2` in the notebook indicates that it excludes the AE sensors affected by the reflected P waves. `Case 1` was without the exclusion, which was obsolated in this study.

7. `p07_solve_lsq_gain/solve_lsqprob_gain_case2.m` and `p07_solve_lsq_gain/solve_lsqprob_direct_and_mixedmodel_case2.m` performs the least square method to obtain the model parameters fitting the P wave amplitude between the synthetic and observed waveforms. Note that we conducted the case studies with the different amplitude models to evaluate the performance. The main model $\hat{A}_{ij} = A_{ij}^{\text{model}} S_i T_j \beta(\omega, \theta)$ is in the `Case 4. mixed model` of `solve_lsqprob_direct_and_mixedmodel_case2.m`.

8. `08_plot_surfaceeffect_result_Case2.ipynb` plots the results for the gain factors and the comparison of the waveforms.