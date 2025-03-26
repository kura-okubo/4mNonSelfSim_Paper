#!/usr/bin/env python
# coding: utf-8

# # Dynamic rupture modeling of gouge patch: compare the M0 and STF
# 
# In this notebook, we compute the seismic moment and source time function obtained from the dynamic rupture modeling.
# We dump the processed data into HDF5.
# We evaluate the source parameters by fitting the synthetic STF in the later notebooks.
# 
# 2024.02.22 Kurama Okubo
# 
# - 2025.1.30 Clean up the notebook for the master plot.
# - 2025.2.3 Update for parameter study
# - 2025.3.22 Update for master case

# Save the source parameter of dynamic rupture models to conduct the grid search

# In[1]:


import os
import shutil
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import ScalarMappable
# get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import pandas as pd
from datetime import timedelta
from tqdm import tqdm
import warnings
import time
from datetime import datetime

import pickle

from scipy.optimize import minimize

from scipy.signal import freqz
from scipy import signal

import h5py
import seaborn as sns

from post_dynrup_func import *

# get_ipython().run_line_magic('load_ext', 'autoreload')
# get_ipython().run_line_magic('autoreload', '2')

plt.rcParams["font.family"] = 'Arial'
# plt.rcParams["font.sans-serif"] = "DejaVu Sans, Arial, Helvetica, Lucida Grande, Verdana, Geneva, Lucid, Avant Garde, sans-serif"
plt.rcParams["font.size"] = 12
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["xtick.major.size"] = 4.75
plt.rcParams["xtick.major.width"] = 0.75
plt.rcParams["xtick.minor.size"] = 3
plt.rcParams["xtick.minor.width"] = 0.4
plt.rcParams["xtick.minor.visible"] = True

plt.rcParams["ytick.direction"] = "in"
plt.rcParams["ytick.major.size"] = 4.75
plt.rcParams["ytick.major.width"] = 0.75
plt.rcParams["ytick.minor.size"] = 3
plt.rcParams["ytick.minor.width"] = 0.4
plt.rcParams["ytick.minor.visible"] = True

plt.rcParams["savefig.transparent"] = False

plt.rcParams['axes.linewidth'] = 0.75


# In[2]:


figdir = "../figure/01_M0andSTF"
if not os.path.exists(figdir):
    os.makedirs(figdir)


# In[3]:


datadir = "../data/01_M0andSTF"
if not os.path.exists(datadir):
    os.makedirs(datadir)


# # 1. Compute seismic moment
# 
# The seismic moment is computed as follows:
# 
# $$ M_0(t) = \mu \int_A u(\mathbf{\xi}, t) dA $$
# 
# $$ = \mu dA \sum_i^{N} u_i $$
# 
# where $dA$ is the area of grid in the simulation, which is uniform in the simulaiton.
# 
# # 2. Compute STF
# 
# The seismic moment is computed as follows:
# 
# $$ \dot{M}_0(t) = \mu \int_A \dot{u}(\mathbf{\xi}, t) dA $$
# 
# $$ = \mu dA \sum_i^{N} \dot{u}_i $$
# 
# where $dA$ is the area of grid in the simulation, which is uniform.
# 

# In[4]:


E = 96e9
nu = 0.246 # metagabbro
mu = E/(2*(1+nu))

a_patch = 4.0e-3
rupturetype = "pulse"
pdcscaling= 0.60 #0.54 #0.475 #0.5 #0.55 # 0.65
bgbeta= 0.35 #0.4 #0.3 #0.5
# gammautry = 0.8

nb_x_elements = 1024 #128
nb_z_elements = 1024 #128

sig_n = 6e6

IfBinaryOutput = True

read_dumpedpickle = False # false if you first run the notebook


# In[5]:


# case study parameter casename
casestr = f"a={a_patch*1e3:.2f}_ruptype={rupturetype}_pdcscaling={pdcscaling:.3f}_sn={sig_n/1e6:.1f}MPa_hatfr=0.3_bgbeta={bgbeta:.2f}"

finame=f"../../preprocess_modelsetup/data/gouge_dynamicrupture_modelparam_{casestr}_dev.csv"

# Read model parameters
df_modelparam = pd.read_csv(finame, index_col=0)

# datadir_root = "../../../uguca/build_v42_masterfitmodel/simulations_main_casestudy"

# datadir_root = "/Volumes/4mGouge_WorkHDD/RuptureSimulation/build_hpcv15_neweventset_dt1e-8/simulations_main_casestudy_hpc"
# datadir_root = "/Volumes/4mGouge_WorkHDD/RuptureSimulation/build_hpcv52_paramstudy_v11/simulations_main_casestudy_hpc_dev"
# datadir_root = "/Volumes/4mGouge_WorkHDD/RuptureSimulation/build_hpcv53_paramstudy_master/simulations_main_casestudy_hpc_paramstudymaster"
# datadir_root = "/Volumes/4mGouge_WorkHDD/RuptureSimulation/build_hpcv53_paramstudy_master/simulations_main_casestudy_hpc_paramstudymaster"
datadir_root = "/Volumes/Okuboetal2025_masterHDD/RuptureSimulation/main_casestudy/build_hpcv61_paramstudy_master_v3/simulations_main_casestudy_hpc_paramsearch"

casestr


# In[6]:


df_modelparam.head()


# In[ ]:

# In[7]:


# parameter study
ifParamStudy = True # change the file name for the parameter study

selectids_list = [24, 50, 52, 72, 129]
# selectids_list = [24]
# selectids_list = [24, 50]

# set the test range
step = 0.005
# delsigma_factor_range = [
# # np.arange(0.68, 0.72+step/2, step=step),
# np.arange(0.7, 0.7+step/2, step=step),
# # np.arange(0.485, 0.525+step/2, step=step),
# np.arange(0.47, 0.53+step/2, step=step),
# # np.arange(0.45, 0.49+step/2, step=step),
# np.arange(0.44, 0.48+step/2, step=step),
# # np.arange(0.405, 0.445+step/2, step=step),
# # np.arange(0.385, 0.425+step/2, step=step),
# np.arange(0.41, 0.44+step/2, step=step),
# np.arange(0.39, 0.42+step/2, step=step)]

# for master v3
delsigma_factor_range = [
np.arange(0.68, 0.72+step/2, step=step),
np.arange(0.465, 0.54+step/2, step=step), 
np.arange(0.44, 0.48+step/2, step=step),
np.arange(0.405, 0.445+step/2, step=step),
np.arange(0.385, 0.425+step/2, step=step)]


delsigma_factor_range

# In[9]:
# select model parameter index
# In[10]:

# delsigma_factorstr = dict() 
# for ii, sid in enumerate(selectids_list):
#     delsigma_factorstr[f"{sid}"] = delsigma_factor_range[ii][delsigma_index[f"{sid}"]]
    
# delsigma_factorstr


# In[ ]:

# In[11]:


expr_id = 87
# model_ids = [24, 50, 52, 72, 129]
# model_ids = [24, 52, 129]
# model_ids = [50, 52, 72, 129]
# model_ids = [24]

#--------------------------------------------------------#
# Parameter search: loop on the model ids and the factors
#--------------------------------------------------------#
# data_paramstudy = {
#     "model_id": [model_ids[0]],
#     "M0_obs": df_gougeevent_selected_modelled["M0"].values[0],
#     "Tw_obs": df_gougeevent_selected_modelled["Tw"].values[0],
#     "M0_model": df_dynrup_sourceparam["M0_bestfit"].values[0],
#     "Tw_model": df_dynrup_sourceparam["Tw_bestfit"].values[0]         
# }

df_paramstudy_all = pd.DataFrame(columns=["model_id", "delsigma_factor", "M0_obs", "Tw_obs", "M0_model", "Tw_model"])


for ii, i_model_ids in enumerate(selectids_list):
    model_ids = [i_model_ids]

    for delsigma_factor in delsigma_factor_range[ii]:

        data_all = dict()

        for model_id in model_ids:
        #model_id = model_ids[0]

            df_modelparam_selected = df_modelparam[df_modelparam.index == model_id]

            if ifParamStudy:
                simulation_name = f"fb03-{expr_id:03d}__{df_modelparam_selected.index[0]:04d}_{casestr}_{delsigma_factor:.4f}"
            else:
                simulation_name = f"fb03-{expr_id:03d}__{df_modelparam_selected.index[0]:04d}_{casestr}"

            print(f"Process {simulation_name}")
            
            df_time = pd.read_csv(os.path.join(datadir_root,simulation_name+".time"), header=None, sep=' ', index_col=0)
            df_coord = pd.read_csv(os.path.join(datadir_root,simulation_name+".coords"), header=None, sep=' ', index_col=None)
            NT = len(df_time)

            # set coordinate
            xcoord = df_coord.loc[:,0].values
            zcoord = df_coord.loc[:,2].values

            x_length = xcoord.max()
            z_length = zcoord.max()

            dgrid = zcoord[1]-zcoord[0]
            dA = dgrid * dgrid
            print(f"Grid size: {dgrid*1e3}[mm]") 

            # read displacement
            read_parameter = "top_disp_0" # consider only the x direction.
            
            
            if IfBinaryOutput:
                D = np.fromfile(os.path.join(datadir_root,simulation_name+f"-DataFiles/{read_parameter}.out"), dtype="float32")
                df_data = pd.DataFrame(data=D.reshape((NT, -1)))
            else:
                df_data = pd.read_csv(os.path.join(datadir_root,simulation_name+f"-DataFiles/{read_parameter}.out"), header=None, sep=' ')        

            # trim the domain
            x_maxwidth=15e-3 # half width of the area to compute the moment; 15mm is large enough to encompass all the slip region
            z_maxwidth=15e-3
            M0_internal_inds_rec = np.where((np.abs(xcoord - dgrid/2 - x_length/2) <= x_maxwidth) & (np.abs(zcoord - dgrid/2 - z_length/2) <= z_maxwidth))

            # extract within the patch    
            rcoord = np.linalg.norm(np.vstack([xcoord - dgrid/2 - x_length/2, zcoord - dgrid/2 - z_length/2]), axis=0)
            M0_internal_inds_patch = np.where(rcoord <= a_patch)
            # M0_internal_inds

        #     h = plt.scatter(xcoord - x_length/2, zcoord - z_length/2, c=rcoord, cmap='viridis')
        #     plt.colorbar(h)

            M0_rec = np.zeros(NT)
            M0_patch = np.zeros(NT)

            for i in range(NT):
                M0_rec[i] = mu * dA * df_data.loc[i, M0_internal_inds_rec].sum()*2 # convert from top displacement to slip by multiplying 2
                M0_patch[i] = mu * dA * df_data.loc[i, M0_internal_inds_patch].sum()*2 # convert from top displacement to slip by multiplying 2

            2    # STF_patch = np.gradient(M0_patch, df_time.values.squeeze())

            # Read velocity
            read_parameter = "top_velo_0" # select the parameter to read
            
            if IfBinaryOutput:
                D = np.fromfile(os.path.join(datadir_root,simulation_name+f"-DataFiles/{read_parameter}.out"), dtype="float32")
                df_data = pd.DataFrame(data=D.reshape((NT, -1)))
            else:
                df_data = pd.read_csv(os.path.join(datadir_root,simulation_name+f"-DataFiles/{read_parameter}.out"), header=None, sep=' ')
                
            STF_rec = np.zeros(NT)
            STF_patch = np.zeros(NT)


            for i in range(NT):
                STF_rec[i] = mu * dA * df_data.loc[i, M0_internal_inds_rec].sum()*2 # convert from top half velocity to the slip velocity by multiplying 2
                STF_patch[i] = mu * dA * df_data.loc[i, M0_internal_inds_patch].sum()*2 # convert from top half velocity to the slip velocity by multiplying 2


            # store the data

            key_t = f"t_{simulation_name}"
            key_M0 = f"M0_rec_{simulation_name}"
            key_STF = f"STF_rec_{simulation_name}"
            key_STF_patch = f"STF_patch_{simulation_name}"


            data_all[key_t] = df_time[1].values
            data_all[key_M0] = M0_rec
            data_all[key_STF] = STF_rec
            data_all[key_STF_patch] = STF_patch
                  
            print(f"{M0_rec[-1]:.3f}J, {M0_patch[-1]:3f}J, {(M0_rec[-1] - M0_patch[-1])/M0_rec[-1]}")
        #     # Dump data
        #     with open(datadir+f'/M0andSTF_all_{tau_r:.2f}.pickle', 'wb') as fo:
        #         pickle.dump(data_all, fo, protocol=pickle.HIGHEST_PROTOCOL)


        # In[12]:


        # time step of numerical simulation

        dt_dynrup = (df_time.values[1]-df_time.values[0])[0]
        print(f"time step of dynamic rupture model is {dt_dynrup*1e6:.4f} μs, {1/dt_dynrup/1e6:.2f} MHz.")


        # In[ ]:





        # # Plot time history of M0
        # 
        # Here we plot the time history of $M_0$ to check the growth of slip.

        # In[ ]:





        # In[13]:


        fig, ax = plt.subplots(1, 1, figsize=(8, 6))

        for model_id in model_ids:
            
            df_modelparam_selected = df_modelparam[df_modelparam.index == model_id]

            if ifParamStudy:
                simulation_name = f"fb03-{expr_id:03d}__{df_modelparam_selected.index[0]:04d}_{casestr}_{delsigma_factor:.4f}"
            else:
                simulation_name = f"fb03-{expr_id:03d}__{df_modelparam_selected.index[0]:04d}_{casestr}"

            key_M0 = f"M0_rec_{simulation_name}"
            
            ax.plot(df_time*1e6, data_all[key_M0], "k-")
            #     ax.plot(df_time*1e6, M0_patch, "r-", label="Dynamic rupture model patch domain")
            # ax.plot(df_time*1e6, cosine_stf, "gray", ls="--",  label=f"Cosine STF: TR={TR*1e6:.1f}μs")

            ax.axhline(df_modelparam_selected.M0_mean.values[0], ls="--", c="k")

        ax.set_xlim([0, 10])
        # ax.set_ylim([-0.02, 1.1])

        ax.set_xlabel("Time [μs]")

        ylabelstr = r"$\mathrm{{M}_0}$"
        ax.set_ylabel("{} [Nm]".format(ylabelstr))

        # ax.legend(loc=0)

        plt.tight_layout()
        # plt.savefig(figdir +f"/M0_all_{casestr}.png", dpi=300, bbox_inches="tight")

        # plt.close()
        # plt.clf()


        # In[ ]:





        # # Plot source time function
        # 
        # We deprecated the comparison to the stacked observation of STF. We compare the modeled STF with the best-fit averaged cosine STF.

        # In[14]:


        # load color dictionary consistent to the plot of the repeated waveforms
        repeated_sensor_lcdict = "OL08" # the color dict is same for all the sensor although separately saved.
        gougepatch_id = "G3"
        with open(f'../../../../ComputeScaling/data/01_plot_gougeevents/lc_dict_{gougepatch_id}_{repeated_sensor_lcdict}.pkl', 'rb') as fi:
            lc_dict = pickle.load(fi)


        # In[15]:


        fig, ax = plt.subplots(1, 1, figsize=(8, 6))

        tvec_model = df_time.values.squeeze()

        i = 0

        maxloc_center = 1.3e-6

        # tvec_data = np.array(fo_stack[f"param/tvec_upsampled_trimmed"])
        # dt_upsampled = tvec_data[1] - tvec_data[0]

        for model_id in model_ids:

            # plot observation data and synthetic STF
            datacase = f"fb03-087__{model_id:04d}"

            # plot synthetic
            df_modelparam_selected = df_modelparam[df_modelparam.index == model_id]
            M0_mean = df_modelparam_selected["M0_mean"].values[0]
            Tw_mean = df_modelparam_selected["Tw_mean"].values[0]
               
            # compute synthetic STF
            tvec_syn = np.linspace(0, Tw_mean, int(Tw_mean/dt_dynrup))
            STF_syn = stf_cosine(tvec_syn, Tw_mean, M0_mean)
            # STF_syn = stf_kupper(tvec_syn, Tw_mean, M0_mean)
            
            ax.plot((tvec_syn+0)*1e6, STF_syn/1e6, ls=":", c=lc_dict[datacase], lw=2.0) # before alignment
            
            # Plot dynamic rupture model
            df_modelparam_selected = df_modelparam[df_modelparam.index == model_id]
            
            if ifParamStudy:
                simulation_name = f"fb03-{expr_id:03d}__{df_modelparam_selected.index[0]:04d}_{casestr}_{delsigma_factor:.4f}"
            else:
                simulation_name = f"fb03-{expr_id:03d}__{df_modelparam_selected.index[0]:04d}_{casestr}"

            key_STF = f"STF_rec_{simulation_name}"
            STF_rec = data_all[key_STF]
            STF_maxarg = np.argmax(STF_rec)
            dt_model = tvec_model[1] - tvec_model[0]
            tshift = (tvec_model[STF_maxarg] - maxloc_center)
            ax.plot((tvec_model - tshift)*1e6, STF_rec/1e6, "-", c=lc_dict[datacase], lw=2, zorder=3)

            # Apply band-pass filter to mimic the observation
            freqmin = 0.1e6
            freqmax = 1e6
            butterworth_order = 3
            filtered_yshift = [0.02, 0.08, 0.12, 0.18 ,0.25]
            b, a = signal.butter(butterworth_order, (freqmin, freqmax), 'bandpass', fs=1/dt_model, output='ba')
            STF_rec_filtered = signal.filtfilt(b, a, STF_rec, method='gust') # using two-way filter Gustafsson’s method
            # ax.plot((tvec_model - tshift)*1e6, STF_rec_filtered/1e6+ filtered_yshift[i], "--", c=lc_dict[datacase], lw=2, zorder=3)

            
            i+=1


        ylabelstr = r"$\dot{M}_0(t)$"
        ax.set_xlabel("Time [μs]")
        ax.set_ylabel("{} [MNm/s]".format(ylabelstr))

        ax.set_xlim([-2, 4.5])
        ax.set_ylim([-0.05, 1.0])
        # ax.set_ylim([-0.05, 0.1])


        # In[ ]:





        # # Compute the source parameters $M_0$ and $T_w$
        # 
        # Here we compute the $M_0$ and $T_w$ by fitting the synthetic STF same as the main analysis of the observations. For a fair comparison, we apply the attenuation factor to mimic the path effect, then apply the low-pass filter similar with the main analysis, and deconvolve the attenuation factor with the water-level. We fit the processed STF with the synthetic cosine STF to estimate the source parameters. Note that we ignore the directivity effect, which modifies the shape of STF in the observations. See the `Others/Synthetictest_STFestimation` for the synthetic test of this process flow.
        # 
        # > We apply the low-pass filter instead of the band-pass filter to avoid the acausal artifacts on the STF. In the main analysis, we applied the band-pass filter on the velocity waveform, and detrend using the polynomials, which mitigated the acausal bump in the STF. The purpose of the band-pass filter is to remove the low-frequency noise on the AE waveforms. Here, since the dynamic rupture model does not show the noise, we selected the low-pass filter at 1MHz same as the upper bound of main analysis to estimate the source parameters.

        # <img src="01_dynrupSTFfit_schematic.png" alt="01_dynrupSTFfit_schematic.png" style="width: 1000px;"/>
        # 

        # 

        # ## Read Q model
        # 

        # In[16]:


        gougepatch_id = "G3"
        Qinv_quart = 50

        df_Qinv_quantile = pd.read_csv("../../../../Calibration/Attenuation/data/df_Qinv_quantile.csv", index_col=0)
        df_Qinv_quantile.head()



        # ## Process the dynamic rupture STF

        # In[17]:


        # Processing parameters
        # upsample the data
        dt_upsampled = 1e-9
        
        zerowin_pre = 10e-6
        zerowin_post = 10e-6
        vp = 6200 #[m/s]
        k_waterlevel = 0.3 # used in the observation analysis

        # Parameters for filtering
        # We use the same filter as the previously analyzed gouge events
        freqmin = 0.1e6 #
        freqmax = 1e6 # 

        butterworth_order = 3


        # fitting STF parameters
        LBA_buffer_winlen = 1.5e-6 # standard value for the LBA buffer

        Tshift_init = 0.0
        pwin_pre = zerowin_pre

        bounds = [(0, 10), (0.1e-6, 10.0e-6), (-1e-6, 1e-6)]
        stf_type = "cosine" #"kupper"
        xatol = 1e-8 #1e-8 #2 
        fatol = 1e-8 #1e-8 #2

        residu_win = [0.5, 0.25]


        # In[18]:


        # Read the source distance datasheet to compute the mean source distance of the event
        trimP_coef_columns = ["rdist", "incidentangle", "dip", "azimuth", "k_M0uz", "TR", "beta_coef_p"]
        df_trimP_coef = pd.DataFrame(columns=trimP_coef_columns)

        AEsensor_list = ["OL23", "OL07", "OL08", "OL22"] # update: we use 4 close sensors

        for stnm in AEsensor_list:
            df_trimP_sensor = pd.read_csv(f"../../../../ComputeScaling/data/02_trim_pwave/trimP_coefficients_{gougepatch_id}_{stnm}.csv", index_col=0)
            df_trimP_coef = pd.concat([df_trimP_coef if not df_trimP_coef.empty else None, df_trimP_sensor])


        # In[ ]:





        # In[19]:


        # # for model_id in model_ids:
        # model_id = model_ids[3]
        # print(f"prcess {model_id}")


        # In[20]:


        lc_debug = sns.color_palette("colorblind")
        lc_debug


        # In[21]:


        df_dynrup_sourceparam = pd.DataFrame(columns=["event_id", "M0_rec", "M0_bestfit", "Tw_bestfit", "Tshift_bestfit"])

        for model_id in model_ids:
        # model_id = model_ids[3]
            print(f"prcess {model_id}")
            
            #------------------------#
            # 1. Zeropadding the STF
            #------------------------#
            
            df_modelparam_selected = df_modelparam[df_modelparam.index == model_id]
            
            if ifParamStudy:
                simulation_name = f"fb03-{expr_id:03d}__{df_modelparam_selected.index[0]:04d}_{casestr}_{delsigma_factor:.4f}"
            else:
                simulation_name = f"fb03-{expr_id:03d}__{df_modelparam_selected.index[0]:04d}_{casestr}"
            
            # simulation_name = f"fb03-{expr_id:03d}__{df_modelparam_selected.index[0]:04d}_{casestr}"
            key_M0 = f"M0_rec_{simulation_name}"
            M0_rec = data_all[key_M0]
            
            key_STF = f"STF_rec_{simulation_name}"
            STF_rec = data_all[key_STF]
            tvec_dynrup = df_time.values.squeeze()
            dt_dynrup = tvec_model[1] - tvec_model[0]
            # plt.plot(tvec_model*1e6, STF_rec/1e6, "-", c=lc_dict[datacase], lw=1, zorder=3)
            
            # upsample the data
            tvec_upsampled = np.arange(tvec_dynrup[0], tvec_dynrup[-1], step=dt_upsampled)
            STF_rec_upsampled = np.interp(tvec_upsampled, tvec_dynrup, STF_rec)
            
            # tpre = -np.arange(dt_dynrup, zerowin_pre, step=dt_dynrup)[::-1]
            # tpost = np.arange(tvec_dynrup[-1]+dt_dynrup, tvec_dynrup[-1]+zerowin_post, step=dt_dynrup)
            tpre = -np.arange(dt_upsampled, zerowin_pre, step=dt_upsampled)[::-1]
            tpost = np.arange(tvec_dynrup[-1]+dt_upsampled, tvec_dynrup[-1]+zerowin_post, step=dt_upsampled)
            tvec_dynrup_padded = np.hstack([tpre, tvec_upsampled, tpost])
            
            post_add = 0
            if np.mod(len(tvec_dynrup_padded), 2) == 1:
                # make tvec length as even
                tvec_dynrup_padded = np.hstack([tvec_dynrup_padded, tvec_dynrup_padded[-1]+dt_upsampled])
                post_add = 1
            
            STF_rec_padded = np.hstack([np.zeros(len(tpre)), STF_rec_upsampled, np.zeros(len(tpost)+post_add)])
            # plt.plot(tvec_dynrup_padded, STF_rec_padded)
            
            #------------------------#
            # 2. Compute attenuation factor
            #------------------------#
            # We use the mean source distance of the four AE sensors to compute the p wave arrival time
            sourcedist_average = df_trimP_coef.loc[f"fb03-{expr_id:03d}__{model_id:04d}"]["rdist"].mean()
            tt_average = sourcedist_average/vp
            print(f"Averaged travel time: {tt_average*1e6:.2f}μs")
            
            Ndata_FFT = len(tvec_dynrup_padded)
            # NFFT = 2**(Ndata_FFT-1).bit_length()
            NFFT = Ndata_FFT # same length of the data for the sake of simplicity
            
            # print(Ndata_FFT, NFFT)
            F_freq = np.fft.rfftfreq(NFFT, d=dt_upsampled)
            Qinv_interp = get_Qinv(F_freq, df_Qinv_quantile.freq.values*1e6, df_Qinv_quantile[f"Qinv_{Qinv_quart}"].values).astype(float)
            
            # plt.loglog(F_freq/1e6, Qinv_interp.real, ".-")
            # plt.xlim([0.1, 2])
            # plt.ylim([0.002, 0.2])
            # plt.xlabel("Frequency [MHz]")
            # plt.ylabel("Attenuation, $Q^{-1}$")
            
            Bomega_interp = np.exp(-np.pi * F_freq * tt_average * Qinv_interp)
            Bomega_wlv = np.maximum(np.abs(Bomega_interp), (k_waterlevel*np.abs(Bomega_interp).max()))
            
            # fig, ax = plt.subplots(1, 1, figsize=(7, 6))
            # ax.loglog(F_freq/1e6, Bomega_interp)
            # ax.loglog(F_freq/1e6, Bomega_wlv)
            # ax.set_xlim([0.06, 2])
            # ax.set_ylim([0.06, 1.2])
            # ax.set_xlabel("Frequency [MHz]")
            # ax.set_ylabel("$B(\omega)$")
            
            #------------------------#
            # 3. Convolve the attenuation factor to mimic the path effect
            #------------------------#
            # compute source spectrum
            F_STF1 = np.fft.rfft(STF_rec_padded, n=NFFT)
            
            #1. convolve the attenuation factor
            STF_Qconvolved = np.fft.irfft(F_STF1 * Bomega_interp).real
            
            #------------------------#
            # 4. Apply the low-pass filter
            #------------------------#
            # b, a = signal.butter(butterworth_order, (freqmin, freqmax), 'bandpass', fs=(1/dt_upsampled), output='ba') # not apply band-pass to mitigate acausal signal
            b, a = signal.butter(butterworth_order, freqmax, 'lowpass', fs=(1/dt_upsampled), output='ba') #
            STF_Qconvolved_filtered = signal.filtfilt(b, a, STF_Qconvolved, method='gust')
            
            
            #------------------------#
            # 5. Deconvolve Bomega with water level
            #------------------------#
            F_STF2 = np.fft.rfft(STF_Qconvolved_filtered, n=NFFT)
            STF_Q_deconvolved = np.fft.irfft(F_STF2/Bomega_wlv).real # divide the spectra by the attenuation factor
            
            #------------------------#
            # 6. Fit the cosine STF to estimate the source parameters
            #------------------------#
            
            # process flow:
            # 1. remove offset at the LBA
            # 2. compute half maximum pulse width to estimate Tw_init
            # 3. search the best-fit source parameters
              
            # pick the LBA as it is
            STF_grad = np.gradient(STF_Q_deconvolved)
                    
            # https://stackoverflow.com/a/3843124
            zero_crossings = np.where(np.diff(np.sign(STF_grad)) > 0)[0]
            min_list = np.array(zero_crossings)
            
            LBA_ind = min_list[np.where(int((zerowin_pre+LBA_buffer_winlen)/dt_upsampled) - min_list > 0)[0][-1]] # search the first bump of STF;
            LBA_amp = STF_Q_deconvolved[LBA_ind]
            LBA_t = tvec_dynrup_padded[LBA_ind]
            
            # remove the offset
            STF_Q_deconvolved_offsetremoved = STF_Q_deconvolved-LBA_amp
            
            # compute Tw_init as the HMPW
            pmax = np.max(STF_Q_deconvolved_offsetremoved)
            pmax_ind = np.argmax(STF_Q_deconvolved_offsetremoved)
            
            halfamp = pmax/2
            
            # search the half pulse width
            halfamp_list = np.where(np.diff(np.sign(STF_Q_deconvolved_offsetremoved - halfamp)))[0]
            
            # for tiny events, skip if we cannot find the LHA or RHA due to low S/N
            LHA_ind = halfamp_list[np.where(halfamp_list - pmax_ind < 0)[0][-1]]
            RHA_ind = halfamp_list[np.where(halfamp_list - pmax_ind > 0)[0][0]] # set HMPW just below the half-maximum amplitude
            HMPW = tvec_dynrup_padded[RHA_ind] - tvec_dynrup_padded[LHA_ind]
            Tw_init = 2*HMPW
            M0_init = 0.5*np.max(STF_Q_deconvolved_offsetremoved)*Tw_init
            Tshift_init = 0.0

            # NOTE: The STF of dynamic rupture is narrow. Thus, the HMPW is not enough for the T_init. 
            # To better and stably fit the cosine STF to the model, we increase the T_init.
            # Update: to avoid the jump in the residual during the grid search caused by the dependency of the initial value,
            # we iterate to find the global minimum in the fitting.
            Tw_init_list = Tw_init * np.linspace(1.0, 1.1, 3)
            res_fun_all = []
            # Tw_init *= 1.05
        
            for Tw_init_test in Tw_init_list:
                x0 = [M0_init, Tw_init_test, Tshift_init]
                # print(x0)
                
                # fit the synthetic STF to the dynamic rupture STF
                res_test = minimize(compute_res, x0, args=(STF_Q_deconvolved_offsetremoved, dt_upsampled, pwin_pre, residu_win, stf_type, False), 
                               method='Nelder-Mead', bounds=bounds, options={"return_all": False, "xatol":xatol, "fatol":fatol, "maxfev":1000})
                res_fun_all.append(res_test.fun)
            
            # use Tw_init with minimum residual 
            Tw_init_best = Tw_init_list[np.argmin(res_fun_all)]
            x0 = [M0_init, Tw_init_best, Tshift_init]
            res = minimize(compute_res, x0, args=(STF_Q_deconvolved_offsetremoved, dt_upsampled, pwin_pre, residu_win, stf_type, False), 
                        method='Nelder-Mead', bounds=bounds, options={"return_all": False, "xatol":xatol, "fatol":fatol, "maxfev":1000})

            
            # x0 = [M0_init, Tw_init, Tshift_init]
            # # print(x0)
            # # fit the synthetic STF to the dynamic rupture STF
            # res = minimize(compute_res, x0, args=(STF_Q_deconvolved_offsetremoved, dt_upsampled, pwin_pre, residu_win, stf_type, False), 
            #                method='Nelder-Mead', bounds=bounds, options={"return_all": False, "xatol":xatol, "fatol":fatol, "maxfev":1000})
            
            # synthesize the estimated STF
            M0_best, Tw_best, tshift_best = res.x
            # print(Tw_best*1e6)
            tvec_syn = np.linspace(0, Tw_best, int(Tw_best/dt_upsampled))
            STF_syn = stf_cosine(tvec_syn, Tw_best, M0_best)
            
            # Debug plot
            fig, ax = plt.subplots(1, 1, figsize=(7, 5.2))
            ax.plot(tvec_dynrup_padded*1e6, STF_rec_padded/1e6, label="original STF", c="k", zorder=2)
            ax.plot(tvec_dynrup_padded*1e6, STF_Qconvolved/1e6, label="Apply B(ω)", c=lc_debug[0])
            ax.plot(tvec_dynrup_padded*1e6, STF_Qconvolved_filtered/1e6, label=f"Apply lowpass \nat {freqmax/1e6:.1f}MHz", c=lc_debug[1])
            ax.plot(tvec_dynrup_padded*1e6, STF_Q_deconvolved_offsetremoved/1e6, "-", 
                    label=f"Deconv B(ω) \nwith waterlevel at {k_waterlevel:.2f}", c="crimson")
            # ax.xlim([-4, 6])
            ax.plot(LBA_t*1e6, 0, "kv", ms=6, zorder=5)
            ax.plot((tvec_syn+tshift_best)*1e6, STF_syn/1e6, "--", label=f"best-fit {stf_type} STF", c="b")
            # ax.plot(tvec_dynrup_padded[zero_crossings]*1e6, np.zeros(len(tvec_dynrup_padded[zero_crossings])), "ks", ms=4)
            # ax.axhline(0)
            ax.legend(loc=0)
            
            props = dict(boxstyle='square', facecolor='white', alpha=1.0)
            
            annot_txt = '\n'.join((
                r'event {}'.format(model_id),
                r'$M_0$={:.2f} Nm'.format(M0_best),
                r'$T_w$={:.2f} μs'.format(Tw_best*1e6),
                r'source dist.={:.1f}mm'.format(sourcedist_average*1e3)))
            
            ax.text(0.05, 0.9, annot_txt, transform=ax.transAxes, fontsize=11,
                    verticalalignment='top', bbox=props)
            
            # ax.set_xlim(np.array([tvec_dynrup_padded[0], tvec_dynrup_padded[-1]]) * 1e6)
            ax.set_xlim([-3, 8])
            ax.set_xlabel("Time [μs]")
            ylabelstr = r"$\dot{M}_0$"
            ax.set_ylabel("{} [MNm/s]".format(ylabelstr))
            ax.grid(True, c=np.array([230, 230, 230])/255, lw=0.25, zorder=-1, which="major")
            ax.set_axisbelow('True')
            
            fig.tight_layout()
            
            # plt.savefig(figdir + f"/dynrupSTFfit_{gougepatch_id}_event{model_id}.png", dpi=80)
            # plt.savefig(figdir + f"/dynrupSTFfit_{gougepatch_id}_event{model_id}.pdf")
            # 
            
            # save the best-fit parameters
            data = {"event_id":[model_id],
                    "M0_rec":[M0_rec[-1]],
                    "M0_bestfit":[M0_best],
                    "Tw_bestfit":[Tw_best],
                    "Tshift_bestfit":[tshift_best]}
            
            df_param = pd.DataFrame.from_dict(data)
            df_dynrup_sourceparam = pd.concat([df_dynrup_sourceparam if not df_dynrup_sourceparam.empty else None, df_param])
            
            plt.clf()
            plt.close()


        # In[22]:


        Tw_init


        # In[23]:


        tshift_best


        # In[24]:


        df_dynrup_sourceparam


        # In[25]:


        # dump the dynamic rupture source parameters to the datasheet
        # df_dynrup_sourceparam.to_csv(f"../data/dynrup_bestfit_sourceparam_{casestr}.csv", index=False, float_format="%12.8g")


        # In[ ]:





        # ## Plot the scaling

        # In[26]:


        df_gougeevent = pd.read_csv(f"../../../../GougeEventCatalog/data/gougeeventcatalog__fb03-{expr_id:03d}__{gougepatch_id}__Q{Qinv_quart}.csv")


        # In[27]:


        Nvalidsensors_thresh = 4
        df_gougeevent_selected = df_gougeevent[df_gougeevent["Nvalidsensors"] == Nvalidsensors_thresh]


        # In[28]:


        print(len(df_gougeevent_selected))
        df_gougeevent_selected.head()


        # In[29]:


        scatter_mc = sns.color_palette("Set1")
        scatter_mc


        # In[30]:


        df_gougeevent_selected_modelled = df_gougeevent_selected[df_gougeevent_selected["event_id"].isin(model_ids)]


        # In[31]:


        # compute regression
        colnames = ["Method", "Intercept", "Slope", "Angle(degrees)", "P-perm(1-tailed)"]
        df_res = pd.read_csv("../../../../ComputeScaling/data/07_loglinearfit/lmodel2_out_regression.txt", sep=' ',
                             names=colnames, skipinitialspace=True, skiprows=1, header=None)
        df_res = df_res.set_index("Method")
        fit_method="MA"

        M0_reg = np.logspace(-3, 1, 101)
        TR_reg = (10**df_res.loc[fit_method, "Intercept"])*(M0_reg**df_res.loc[fit_method, "Slope"])


        # In[32]:


        fig, ax = plt.subplots(1, 1, figsize=(5.7, 6))

        scatter_kws0 = {"s": 0, "edgecolors": "k", "zorder": 1, "alpha": 0.9}
        scatter_kws1 = {"s": 90, "edgecolors": "k", "zorder": 1, "alpha": 0.9}
        line_kws = {"color": "crimson", "zorder": -3}

        tc = [""]
        mctype = ["o", "d", "s", "v"]

        labelflag = 0

        ifPlotMain = True
        ifPlotModeled = True
        ifPlotRegression = True

        # Compute standard error
        standarderror_factor = np.sqrt(Nvalidsensors_thresh)

        mainmarkersize = 7

        if ifPlotMain:
            ax.errorbar(df_gougeevent_selected["M0"].values, df_gougeevent_selected["Tw"].values*1e6, 
                    yerr = df_gougeevent_selected["Tw_std"].values*1e6/standarderror_factor, xerr = df_gougeevent_selected["M0_std"]/standarderror_factor,
                    capsize=0, fmt='o', markersize=mainmarkersize, color=scatter_mc[0], lw=1, markeredgecolor = "black", label="Mean of four AE sensors", zorder=3)

        if ifPlotModeled:
            ax.errorbar(df_gougeevent_selected_modelled["M0"].values, df_gougeevent_selected_modelled["Tw"].values*1e6, 
                yerr = df_gougeevent_selected_modelled["Tw_std"].values*1e6/standarderror_factor, xerr = df_gougeevent_selected_modelled["M0_std"]/standarderror_factor,
                capsize=0, fmt='^', markersize=mainmarkersize, color=scatter_mc[1], lw=1, markeredgecolor = "black", label="Mean of four AE sensors", zorder=3)

        # plot the best-fit regression
        if ifPlotRegression:
            ax.plot(M0_reg, TR_reg, c=scatter_mc[0], lw=1.0, zorder = 2)


        # Plot scaling of dynamic rupture model
        ax.plot(df_dynrup_sourceparam["M0_bestfit"].values, df_dynrup_sourceparam["Tw_bestfit"].values*1e6, "*", ms=12, mfc=scatter_mc[5], mec="k",
                label="Dynamic rupture model", zorder=4)
                   
                # yerr = df_gougeevent_selected["Tw_std"].values*1e6/standarderror_factor, xerr = df_gougeevent_selected["M0_std"]/standarderror_factor,
                # capsize=0, fmt='o', markersize=mainmarkersize, color=scatter_mc[0], lw=1, markeredgecolor = "black", label="Mean of four AE sensors", zorder=3)



        xlimit_scaling = [0.004, 3] #10] # check 1/3
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(xlimit_scaling)
        ax.set_ylim([1.2, 3.8]) #10]) # update 2025/1/23 # check 1/3

        # ax.set_yticks([1, 2, 3, 4])
        # ax.set_yticklabels([1, 2, 3, 4])
        ax.set_yticks([2, 3, ]) # update 2025/1/23
        ax.set_yticklabels([2, 3, ]) # update 2025/1/23

        ax.set_xlabel(r"$M_0$ [Nm]")
        ax.set_ylabel(r"$T_w$ [μs]")

        ax.grid(True, c=np.array([230, 230, 230])/255, lw=0.25, zorder=-1)
        ax.set_axisbelow('True')

        plt.tight_layout()

        # plt.savefig(figdir+f"/preliminary_dynrupscaling_{casestr}_{nb_x_elements}.png",  dpi=200)
        # plt.savefig(figdir+f"/preliminary_dynrupscaling_{casestr}_{nb_x_elements}.pdf")


        # In[ ]:

        # parameter study: dump the estimated source parameters

        data_paramstudy = {
            "model_id": [model_ids[0]],
            "delsigma_factor": delsigma_factor,
            "M0_obs": df_gougeevent_selected_modelled["M0"].values[0],
            "Tw_obs": df_gougeevent_selected_modelled["Tw"].values[0],
            "M0_model": df_dynrup_sourceparam["M0_bestfit"].values[0],
            "Tw_model": df_dynrup_sourceparam["Tw_bestfit"].values[0]         
        }

        df_paramstudy = pd.DataFrame.from_dict(data_paramstudy)
        df_paramstudy_all = pd.concat([df_paramstudy_all if not df_paramstudy_all.empty else None, df_paramstudy])
            

# dump the parameter study dataframe
df_paramstudy_all.to_csv(f"../data/aux01_dynrup_paramstudy_{casestr}.csv", index=False, float_format="%12.8g")
print(f"parameter study {casestr} successfully done.")

# # In[ ]:
