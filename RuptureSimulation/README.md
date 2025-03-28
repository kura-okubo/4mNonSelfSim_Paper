# Dynamic Rupture Simulation of Gouge Events

This directory contains the dynamic rupture simulation setup and related files.

## Software

We used the Spectral Boundary Integral Method (SBIEM)-based software [UGUCA](https://gitlab.com/uguca/uguca) for dynamic rupture modeling on gouge patches.

We utilized the stable version `1.0.0` from GitLab and extended the friction interface to implement self-healing with the linear slip-weakening law.

The extended version of the software and input files for grid-search and master cases are available in a separate repository: [**4mNonSelfSim_UGUCA**](https://github.com/kura-okubo/4mNonSelfSim_UGUCA).

> **Reference:**
> Kammer, D.S., Albertini, G., & Ke C.-Y. (2021). "UGUCA: A spectral-boundary-integral method for modeling fracture and friction." *SoftwareX*, 15, 100785.

## Environment

The simulations were developed using a Mac mini M1 (16GB RAM, macOS v15.3.2) and the HPC cluster at NIED.

### Package Versions

| Package   | HPC Version |
|-----------|------------|
| `gcc`     | 8.1.0      |
| `cmake`   | 3.29.0-rc1 |
| `OpenMPI` | 4.0.2a1    |
| `FFTW`    | 3.3.10     |

## Computing the Kernel

We precomputed the kernel with a Poisson's ratio of $\nu = 0.246$. Due to file naming conventions, the kernel files were renamed as `nu0.24_hXX.txt`.

```
nu      = 0.246
pstress = False
dt      = 0.05
tcut    = 100
```

> **Note for Apple Silicon Mac Users:** See [this link](https://gitlab.com/uguca/uguca/-/merge_requests/29) for instructions on computing the kernel using `multiprocessing`.

## Contents

### Preprocessing

[**preprocess_modelsetup**](./main_casestudy/preprocess_modelsetup)

- [01_gougepatch_dynrup_main_casestudy.ipynb](./main_casestudy/preprocess_modelsetup/code/01_gougepatch_dynrup_main_casestudy.ipynb): Generates the input file for local test simulations.
- [02_gougepatch_dynrup_main_casestudy_hpc_paramgridsearch.py](./main_casestudy/preprocess_modelsetup/code/02_gougepatch_dynrup_main_casestudy_hpc_paramgridsearch.py): Generates input files for grid search of initial stress levels.
- [03_gougepatch_dynrup_main_casestudy_hpc_master.ipynb](./main_casestudy/preprocess_modelsetup/code/03_gougepatch_dynrup_main_casestudy_hpc_master.ipynb): Generates input files for the best parameter set, including cases with and without self-healing friction. Also includes a few-time-step model to visualize initial parameters.
- [04_generate_latex_initparam_table.ipynb](./main_casestudy/preprocess_modelsetup/code/04_generate_latex_initparam_table.ipynb): Outputs a LaTeX table of the best-fit dynamic rupture model parameters.

> The input files used for the main analysis are located in [**4mNonSelfSim_UGUCA**](https://github.com/kura-okubo/4mNonSelfSim_UGUCA).

### Postprocessing

[**postprocess_dynrup**](./main_casestudy/postprocess_dynrup)

- [01_compute_M0andSTF.ipynb](./main_casestudy/postprocess_dynrup/01_compute_M0andSTF.ipynb): Initial check of the simulation results. Computes the STF time history, evaluates source parameters by fitting a cosine STF, and exports the data.
- [02_plot_initialcondition.ipynb](./main_casestudy/postprocess_dynrup/02_plot_initialcondition.ipynb): Plots the initial stress state based on simulation output to verify correct implementation.
- [03_plot_snapshots_mastercase.ipynb](./main_casestudy/postprocess_dynrup/03_plot_snapshots_mastercase.ipynb): Plots snapshots of the dynamic rupture simulation, including traction history at a point. Requires both self-healing and non-self-healing cases for comparison.
- [04_plot_sheartractionhistory.ipynb](./main_casestudy/postprocess_dynrup/04_plot_sheartractionhistory.ipynb): Plots the shear traction history.
- [05_plot_masterSTF.ipynb](./main_casestudy/postprocess_dynrup/05_plot_masterSTF.ipynb): Plots the source time function of the dynamic rupture model.

### Additional Auxiliary Scripts

- [aux01_compute_M0andSTF_paramstudy.py](./main_casestudy/postprocess_dynrup/code/aux01_compute_M0andSTF_paramstudy.py): Processes simulation results to evaluate source parameters, aiding in identifying the best-fit model.
- [aux01_paramstudy_searchbestfitparam.ipynb](./main_casestudy/postprocess_dynrup/code/aux01_paramstudy_searchbestfitparam.ipynb): Visualizes grid search results.
- [aux03_plot_crosssection_slip.ipynb](./main_casestudy/postprocess_dynrup/code/aux03_plot_crosssection_slip.ipynb): Plots the cross-section profile of slip and shear traction.

