#!/usr/bin/env python
# coding: utf-8

# # Dynamic rupture modeling of gouge patch: Case study
# 
# In this notebook we configure the initial conditions of the dynamic rupture modeling associated with the gouge-mediated seismic event.
# 
# 2024.02.12 Kurama Okubo
# 
# - 2024.02.19 update for linea_coulomb_friction_law
# - 2024.02.20 update for the contrast between $\hat{\tau}_{patch}$ and $\hat{\tau}_{background}$ to arrest the rupture.
# - 2024.02.21 append all the parameter in the input file
# - 2024.02.22 update for case study 1. $\tau_r=0$, 2. $\tau_r=0.6$MPa.
# - 2024.03.07 update for the case of patch expansion
# - 2024.03.15 update for master case study
# - 2024.03.18 update to set dynamic_excess for the case of tau_r=0.6MPa
# - 2024.03.27 update the patch margin from 0.1 mm to 0.08 mm to harmonize with the grid size of 0.04 mm.
# - 2024.05.06 update for advanced rupture model: generate the input files for smooth nucleation model and the rapid + stress free model.
# - 2024.05.15 update for the nucleation by decreasing fp. It did not work, so deprecated.
# 
# updated v2: master casestudy
# - 2024.06.11 update to implement master casestudy for the rupture type (crack-like or self-healing pulse-like) and the scaling exponent of Dc.
# - 2024.09.03 update for the new set of non-self-similar events obtained by the stacking of multiple AE sensor. We also removed the dependency of the Energy budget precalculation.
# 
# updated v3: master casestudy with the merged catalog
# - 2024.12.18 update for the new gouge event catalog.
# - 2025.1.23 update the model with $\sigma_n$=6MPa.
# - 2025.1.29 update for the master casestudy
# - 2025.1.31 update for init parameter input file
# - 2025.2.2 update running parameter study of the initial stress condition.
# - 2025.3.22 update for master plot

# # Summary of this notebook
# 
# This notebook generates the input file for the dynamic rupture simulation with uguca.
# 
# The parameters to search in the case study is $a_{patch}$, rupture type, and $p$ as the trial scaling exponent of Dc.
# For the single run of this notebook produces a set of input files with different events with the predefined parameters above.
# 
# To produce different cases, re-run the notebook with setting different parameters.

# In[319]:

import os
import sys
import shutil
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
# get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import pandas as pd
from datetime import timedelta
from tqdm import tqdm
import warnings
import time
from datetime import datetime

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

plt.rcParams["savefig.transparent"] = True

plt.rcParams['axes.linewidth'] = 0.75


# In[320]:


figdir = "../figure"
if not os.path.exists(figdir):
    os.makedirs(figdir)


# In[321]:


datadir = "../data"
if not os.path.exists(datadir):
    os.makedirs(datadir)


# # Set case study variables

# In[322]:


a_patch = 4.0e-3 # patch radius without margine
rupturetype = "pulse" # "crack": without self-healing or "pulse": with self-healing
p_dcscaleexp = 0.60 #0.55 #0.54 #0.55 #0.54 #0.56 #0.555  # #0.475 #0.55 #0.575 #0.65 #0.8 #0.65 #0.7

IfInitparam = False # True to make the input files with the short time duration to output the initial condition


# In[323]:


casestudy_name = f"a={a_patch*1e3:.2f}_ruptype={rupturetype:s}_pdcscaling={p_dcscaleexp:.3f}"
print(casestudy_name)


# # Theory to derive the model parameters
# 
# The key parameters are the stress drop and $Dc$, which are unknown due to the causal loop: the "true" stress drop and $D_c$ can be estimated from fitting the data to the dynamic rupture simulation, whereas we need first set them to run the simulation. Ideally, the iterative approch can solve this loop. However, it is way too much in our purpose as we already have apriori information on the source paramters. To solve this loop, we set the trial parameters inferred from the gouge patch size and the seismic moment. Note that we set the gouge patch size as the free parameter to answer the question if this affects the conclusions.
# 
# We set the stress drop as followings:
# 
# $$ \bar{ \Delta \sigma } _{\text{try}}^{(i)} = \dfrac{7}{16} \dfrac{M_0^{(i)}}{a^3_{\text{patch}}}, $$
# 
# where $(i)$ indicates the parameter of the $i$th event. We use $a_{\text{patch}} = 4.0$ mm + $0.08$ mm of the patch margin for the initial condition.
# 
# $$ \bar{u}_{\text{try}}^{(i)} = \dfrac{M_0^{(i)}}{\mu \pi a^2_{\text{patch}}} $$
# 
# where $\bar{u}_{\text{try}}^{(i)}$ is the trial averaged slip.
# 
# We set the initial shear stress as the fraction of peak friction:
# 
# $$ \tau_0 = cf_p\sigma_n, $$
# 
# where $c$ determines the coefficient of initial traction fraction. Then, we set the $f_p$ as follows:
# 
# $$ \Delta\sigma = cf_p\sigma_n - f_r\sigma_n $$
# $$ f_p = \dfrac{1}{c}\left[ \dfrac{s\Delta\sigma}{\sigma_n} + f_r \right] $$ 
# 
# Here, $\Delta\sigma$ is controlled with the `dynamic_excess`, $s$, as $ \Delta\sigma = s \bar{ \Delta \sigma } _{\text{try}}^{(i)} $ to control the amplitude of STF.
# 
# To estimate the $D_c$ as follows:
# 
# $$ D_{c}^{(i)} =  \dfrac{D_c^{\text{min}}}{\min \left\{ \bar{u}_{\text{try}} \right\}^{p} } \bar{u}_{\text{try}}^{(i)p},$$
# 
# $D_c^{\text{min}}$ indicates the best-fit $D_c$ for the case of minimum gouge event, which we find by trial and error. $p$ is the scaling exponent with the trial slip.
# 
# When the trial parameters are validated with the observations, i.e., source time functions with non-self-similarity, we can evaluate the "true" paramters such as $\bar{ \Delta \sigma } _{\text{true}}^{(i)}$ and $\bar{u}_{\text{true}}^{(i)}$ as well as the true source region $A_{\text{true}}$.
# 
# **NOTE:**
# 
# The way to define $D_c$ has been updated from our previous work. We used the energy budget such as Dc = (1-$\eta_R$)$\bar{u}$, but it has an assumption of the radiated energy inferred from the kinematic source model. As the direct observation of $E_R$ is unstable and difficult, we changed the derivation as described above. In this metric, we consider the increase in the  $\bar{u}_{\text{try}}^{(i)p}$ from the $D_c$ to model the minimum event, which we predefine with a emperical constant value. 

# # Set model parameters

# In[324]:


# Elastic constant
E = 96e9
rho = 2980
nu = 0.246 # metagabbro
mu = E/(2*(1+nu))

print(f"E, mu, rho = {E:.4g} {mu:.4g} {rho:.4g}")

R_patch = a_patch #4e-3 # gouge patch radius 
R_margin = a_patch+0.08e-3 #4.08e-3 #4.1e-3 #5e-3 # This is the outer bound of the stress margin, used for the input parameter of simulation
# R_margin = a_patch

# Set the range of average GIIC
A_patch = np.pi * R_patch**2

hat_sn_patch = 6e6 #8e6 # normal stress on gouge patch
hat_sn_background = 2e6 # normal stress on background region

hat_fr_patch = 0.3 # fixed the residual friction level as an assumption
hat_tau_r_patch = hat_sn_patch * hat_fr_patch


# ## Recompute case study parameters

# In[325]:


gougepatch_id = "G3" # to set output filename
denoise_method = "detrend"
Qinv_quart = 50
k_waterlevel = 0.3
expr_id = 87

foname_mean = f"../../../../ComputeScaling/data/05_STFstats/SourceParam_meanstd_fb03-{expr_id:03d}_{gougepatch_id}_wlv_{k_waterlevel:.2f}_denoisemethod_{denoise_method.lower()}.csv"
df_stats = pd.read_csv(foname_mean, index_col=2)


# In[326]:


df_stats.head()


# In[327]:


Qinv_quart = "50"
NvalidSensor = 4

df_stats_selected = df_stats[(df_stats["Qinv_quart"] == Qinv_quart) & (df_stats["Nvalidsensors"] >= NvalidSensor)].copy()
df_stats_selected.head()
print(f"Num. event = {len(df_stats_selected)}")


# In[328]:


df_stats_selected["Tw_mean"] * 1e6


# In[329]:


def M02Mw(M0):
    """
    convert from M0 to Mw
    """
    return (np.log10(M0) - 9.105) * 2.0 / 3.0


# In[330]:


df_stats_selected.loc[:, "Mw_mean"] = df_stats_selected.apply(lambda x: M02Mw(x.M0_mean), axis=1)


# In[331]:


M0_stats = (df_stats_selected['M0_mean'].mean(), df_stats_selected['M0_mean'].std(), df_stats_selected['M0_mean'].min(), df_stats_selected['M0_mean'].max())
Mw_stats = (df_stats_selected['Mw_mean'].mean(), df_stats_selected['Mw_mean'].std(), df_stats_selected['Mw_mean'].min(), df_stats_selected['Mw_mean'].max())


# In[332]:


print("M0 stats: {0[0]:.3f} ± {0[1]:.3f} Nm with the range from {0[2]:.3f} to {0[3]:.3f} Nm.".format(M0_stats))
print("Mw stats: {0[0]:.3f} ± {0[1]:.3f} with the range from {0[2]:.3f} to {0[3]:.3f}.".format(Mw_stats))


# In[333]:


# make a dataframe for the dynamic rupture parameters
df_dynparam = df_stats_selected[["M0_mean", "Tw_mean", "Mw_mean"]].copy()
df_dynparam.head()


# In[334]:


df_dynparam.loc[:, "hat_sn_patch"] = hat_sn_patch
df_dynparam.loc[:, "hat_sn_background"] = hat_sn_background


# In[335]:


#1. compute trial stress drop
df_dynparam["delsig_withmargin_try"] = df_dynparam.apply(lambda x: (7/16) * (x.M0_mean/R_margin**3), axis=1)

#2. compute trial average slip
df_dynparam["slip_try"] = df_dynparam.apply(lambda x: x.M0_mean/(mu * np.pi * R_margin**2), axis=1)


# In[336]:


df_dynparam.sort_values("Mw_mean")


# ## Set the frictional parameters

# |  |  $\sigma_n$  | $\tau_0$ | $f_p$ | $f_r$ |
# | ---- | ---- | ---- | ---- | ---- | 
# | nucleation zone | $\alpha \hat{\sigma}_n$ | Gaussian distribution | $f_p$ | $\hat{f}_r$ |
# | gouge patch zone | $\hat{\sigma}_n$ | 0.925*$\tau_p$ | $f_p$ | $\hat{f}_r$ |
# | stress margin | 0 | 0 | 0 | 0 |
# | background region | $\hat{\sigma}_n^{background}$ | $\beta \tau_r$ | $f_p = \hat{f}_r$| $\hat{f}_r$ |
# 
# 

# In[ ]:





# In[337]:


# # Master unirateral + best case

R_nuc = 2.5e-3 #1.5e-3 # nucleation radius 
A_nuc = np.pi * R_nuc**2

nuc_x = -(R_patch - R_nuc) #-(0.5*R_patch) # x coordinate of the center of the nucleation area

# Parameter of the background zone 
print(f"hat_fr_patch={hat_fr_patch}")
hat_fp_background = 0.4 #0.3 # estimated from macroscopic friction value
hat_fr_background = 0.4 #0.3

nuc_normalstress_alpha = 1.0 # amplication factor of the normal stress on the nucleation zone

stressbackground_beta = 0.35 #0.3 #0.4 #0.3 # factor to define the background stress level; this decides the strength of barrier

#--- set the self-healing parameter---#
if rupturetype=="crack":
    hat_ds_factor_rapidnuc_nuc = 10000 # factor of slip-strengthening distance in the nucleation zone
    hat_ds_factor_rapidnuc_patch =  10000 # factor of slip-strengthening distance in the patch area
    
elif rupturetype=="pulse":
    hat_ds_factor_rapidnuc_nuc = 5.5 #6.0 #6.3 #6.5 #5.5  # factor of slip-strengthening distance in the nucleation zone
    hat_ds_factor_rapidnuc_patch = 5.5 #6.0 #6.3 #6.5 #5.5 # factor of slip-strengthening distance in the patch area
    
else:
    raise ValueError(f"rupturetype {rupturetype} not defined.")

hat_ds_factor_rapidnuc_background = 10000 # to avoid the slip-strengthening for the background region

# initialstress_fraction = 0.925 #0.9 #0.9 # initial shear stress is initialstress_fraction*sn*fp

c_nucexcess = 0.02 #0.025 #0.015 # #0.05 #0.05 # the percentage of the excess of the initial shear stress tau0_{nuc}^{max} = (1+c)taup

casename = casestudy_name+"_sn={:.1f}MPa_hatfr={:.1f}_bgbeta={:.2f}".format(hat_sn_patch/1e6, hat_fr_patch, stressbackground_beta)

print(casename)


# In[338]:


# set initial stress fraction on the events
df_dynparam.loc[:, "initialstress_fraction"] = 0.925 #0.875 #0.925


# change the initial stress fraction on the small events
df_dynparam.loc[24, "initialstress_fraction"] = 0.98 #0.94 # large initialstress_fraction 

# df_dynparam.loc["fb03-087__0035", "initialstress_fraction"] = 0.94 #0.925 #0.975


# # We unified the initialstress_fraction with different events
# if rupturetype=="pulse":
#     df_dynparam.loc["fb03-087__0036", "initialstress_fraction"] = 0.942 #0.945 #0.975 # large initialstress_fraction 
#     df_dynparam.loc["fb03-087__0035", "initialstress_fraction"] = 0.94 #0.925 #0.975


# elif rupturetype=="crack":

#     df_dynparam.loc["fb03-087__0036", "initialstress_fraction"] = 0.942 #0.945 #0.975 # large initialstress_fraction 
#     df_dynparam.loc["fb03-087__0035", "initialstress_fraction"] = 0.94 #0.925 #0.975


# In[339]:


df_dynparam.loc[[24, 50, 52, 72, 129], :]


# <!-- ### Compute equivalent delsigma_factor
# 
# We keep the $D_c^{min}$ with `initialstress_fraction` $c$  = 0.925 with `delsigma_factor` $s$ = 0.7. To keep it, we set the delsigma_factor as follows:
# 
# 
# $$ D_c^{min} = \dfrac{2G_{IIC}^{min}}{\sigma_n \left[\dfrac{1}{c}[\dfrac{s\Delta \sigma}{\sigma_n} + f_r] - f_r\right]} $$
# 
# $$s = \dfrac{\sigma_n}{\Delta \sigma} \left[ c [\dfrac{2G_{IIC}^{min}}{D_c^{min}\sigma_n} + f_r] -f_r \right] $$ -->

# In[340]:


# GIICmin = 0.0064519416 # obtained with c=0.925 and s=0.7
# Dcmin = 2.3551006e-08
# c_min = 0.95

# min_gougeid = 128 # gouge event to use as the minimum event in the set of non-self-similar events 

# df_dcmin = df_dynparam[df_dynparam.index == f"fb03-087__{min_gougeid:04d}"]
# dcmin_slip = df_dcmin.slip_try.values[0]

# delsigma_factor_min_fixed = (df_dcmin.hat_sn_patch/df_dcmin.delsig_withmargin_try) * (c_min * (2*GIICmin / (Dcmin*df_dcmin.hat_sn_patch) + hat_fr_patch) - hat_fr_patch)
# print(f"initalstress_fraction={c_min} corresponds to delsigma_factor = {delsigma_factor_min_fixed.values[0]}")


# In[ ]:





# In[ ]:





# ## set dynamic excess

# In[341]:


df_dynparam.loc[:, "delsigma_factor"] = 0.6 # initialize the delsigma


# In[342]:

#Update parameter study 2025.2.2
# loop the delsigma_factor with the event ids
# Note: fix event 24 as minimum event to determine the rest of events 
selectids_list = [24, 50, 52, 72, 129]
# selectids_list = [50]
# selectids_list = [24]

# set the test range
step = 0.005
# delsigma_factor_range = [
# np.arange(0.68, 0.72+step/2, step=step),
# # np.arange(0.695, 0.695+step/2, step=step),
# # np.arange(0.485, 0.525+step/2, step=step),
# np.arange(0.47, 0.53+step/2, step=step), # Next, potentially run from 0.4
# # np.arange(0.45, 0.49+step/2, step=step),
# np.arange(0.44, 0.48+step/2, step=step),
# np.arange(0.405, 0.445+step/2, step=step),
# np.arange(0.385, 0.425+step/2, step=step)]

# 2025.2.26 potential new set of parameter search with a interval of 0.04
# grid search for 5 cases
delsigma_factor_range = [
np.arange(0.68, 0.72+step/2, step=step),
np.arange(0.465, 0.54+step/2, step=step), 
np.arange(0.44, 0.48+step/2, step=step),
np.arange(0.405, 0.445+step/2, step=step),
np.arange(0.385, 0.425+step/2, step=step)]


if IfInitparam:
    inputfileoutdir = f"../../../../../4mNonSelfSim_UGUCA/simulations_main_casestudy_hpc_paramsearch/gouge_rupture_inputfiles_{casestudy_name}_initcondition"
else:
    # inputfileoutdir = f"../../../uguca/simulations_main_casestudy_hpc/gouge_rupture_inputfiles_{casestudy_name}"
    inputfileoutdir = f"../../../../../4mNonSelfSim_UGUCA/simulations_main_casestudy_hpc_paramsearch/gouge_rupture_inputfiles_{casestudy_name}"

# remove the previous case study
if os.path.exists(inputfileoutdir):
    shutil.rmtree(inputfileoutdir)

if not os.path.exists(inputfileoutdir):
    os.makedirs(inputfileoutdir)
    
for ii, selectid in enumerate(selectids_list):
# selectid = selectids_list[0]

    for delsigma_factor_range_event_case in delsigma_factor_range[ii]:
        print(f"process {ii}, factor={delsigma_factor_range_event_case}")
            
        if rupturetype=="pulse":
            df_dynparam.loc[24, "delsigma_factor"] = 0.695 #0.65 # fix the minimum event
            # df_dynparam.loc[50, "delsigma_factor"] = 0.51 #0.515 #0.525 #0.5
            # df_dynparam.loc[52, "delsigma_factor"] = 0.465 #0.475 #0.45
            # df_dynparam.loc[72, "delsigma_factor"] = 0.42 #0.425 #0.425
            # df_dynparam.loc[129, "delsigma_factor"] = 0.4 #0.4 #0.4
            df_dynparam.loc[selectid, "delsigma_factor"] = delsigma_factor_range_event_case #0.4 #0.4

        elif rupturetype=="crack":
            df_dynparam.loc[24, "delsigma_factor"] = 0.695 #0.695 #0.65
            # df_dynparam.loc[50, "delsigma_factor"] = 0.525 #0.5
            # df_dynparam.loc[52, "delsigma_factor"] = 0.475 #0.45
            # df_dynparam.loc[72, "delsigma_factor"] = 0.425
            # df_dynparam.loc[129, "delsigma_factor"] = 0.4
            df_dynparam.loc[selectid, "delsigma_factor"] = delsigma_factor_range_event_case #0.4 #0.4


        # In[ ]:





        # $f_p$ is determined as follows:
        # 
        # $$ f_p = \dfrac{1}{c}\left[ \dfrac{s\Delta\sigma}{\sigma_n} + f_r \right] $$ 
        # 

        # In[343]:


        # compute fp_patch
        # df_dynparam.loc[:, "fp_patch"] = df_dynparam.apply(lambda x: (dynamic_excess*x.delsig_withmargin_try/hat_sn_patch) + hat_fr_patch, axis=1)
        # df_dynparam.loc[:, "fp_patch"] = df_dynparam.apply(lambda x: ((dynamic_excess*x.delsig_withmargin_try/hat_sn_patch) + hat_fr_patch)/initialstress_fraction, axis=1)
        # df_dynparam.loc[:, "fp_patch"] = df_dynparam.apply(lambda x: ((x.delsigma_factor*x.delsig_withmargin_try/hat_sn_patch) + hat_fr_patch)/initialstress_fraction, axis=1) # flexible 
        df_dynparam.loc[:, "fp_patch"] = df_dynparam.apply(lambda x: ((x.delsigma_factor*x.delsig_withmargin_try/hat_sn_patch) + hat_fr_patch)/x.initialstress_fraction, axis=1) # variable initialstress_fraction  

        df_dynparam.loc[:, "hat_fr"] = hat_fr_patch
        df_dynparam.loc[:, "hat_fp_background"] = hat_fp_background
        df_dynparam.loc[:, "hat_fr_background"] = hat_fr_background


        # In[ ]:





        # In[344]:


        df_dynparam.plot.scatter(x="M0_mean", y="fp_patch" )


        # In[345]:


        #3. compute dc

        # compute u_min and dc_min
        min_gougeid = 24 # gouge event to use as the minimum event in the set of non-self-similar events 

        df_dcmin = df_dynparam[df_dynparam.index == min_gougeid]
        dcmin_slip = df_dcmin.slip_try.values[0]
        print(f"min utry = {dcmin_slip*1e6:.3g} μm.")
        # dc_min = 1e-07 # 9.42e-8 we found the best-fit minimum dc by trial and error inferred from kinematic source energy based values;


        # In[346]:


        df_dynparam[df_dynparam.index == min_gougeid]


        # In[ ]:





        # In[347]:


        # Set dc_min from fracture energy scaling
        # load the data
        sei_types={'Reference M0, a':'category'}
        df_seis = pd.read_csv('../../../../Others/EnergyBudget/fracture_energy/merged_data-seismology.csv',sep=';',encoding='cp1252',dtype=sei_types)

        df_m2014 = df_seis[df_seis["Reference M0, a"] == 'McLaskey et al., 2014']
        df_s2019 = df_seis[df_seis["Reference M0, a"] == 'Selvadurai, 2019']
        df_y2014 = df_seis[df_seis["Reference M0, a"] == 'Yoshimitsu et al., 2014']
        df_s2003 = df_seis[df_seis["Reference M0, a"] == 'Sellers et al., 2003']


        # ### Update: the estimation of $G_{IIC}^{syn}$

        # $$ G_{IIC}^{syn} (\delta) = 10^{p_0 + p_1 \log_{10}\delta} $$
        # 
        # We decrease the synthetic slope of $G_{IIC}$ as follows:
        # $$ G_{IIC}^{syn, modified} (\delta) = c 10^{p_0 + p_1 \log_{10} \delta } $$
        # $$ \log_{10} G_{IIC}^{syn, modified} (\delta) = \log_{10} c + p_0 + p_1 \log_{10}\delta $$

        # In[ ]:





        # In[348]:


        from scipy.optimize import curve_fit

        # use the references with the rock sample
        slip_all = []
        GIIC_all = []
        for i, df in enumerate([df_s2003, df_m2014, df_y2014]):
            df = df[df["value (J/m^2)"].astype(float)>0] # select positive GIIC 
            slip_all = np.append(slip_all,  df["value (m)"].astype(float))
            GIIC_all = np.append(GIIC_all,  df["value (J/m^2)"].astype(float))

        # regression with the scaling model
        # ref: https://stackoverflow.com/a/3433503

        popt, pcov = curve_fit(lambda x,a,b: a+b*x ,  np.log10(slip_all),  np.log10(GIIC_all),  p0=(1, 2))

        # UPDATE: decrease the interseption of slope to reproduce the smallest event
        GIIC_slope_intercept_factor = 0.415 #0.425 #0.6 # 0.9

        popt[0] += np.log10(GIIC_slope_intercept_factor)

        slip_syn = np.logspace(-8, -6, 11)

        GIIC_syn = 10**(popt[0] + popt[1]*np.log10(slip_syn))


        # define the function to estimate GIIC from scaling
        def get_GIIC(slp):
            return 10**(popt[0] + popt[1]*np.log10(slp))


        # In[349]:


        popt


        # In[350]:


        df_dcmin


        # In[351]:


        dc_min = 2*get_GIIC(dcmin_slip) / ((df_dcmin.fp_patch-df_dcmin.hat_fr)*df_dcmin.hat_sn_patch)
        print(f"GIIC min: {get_GIIC(dcmin_slip):12.8g}, dc min: {dc_min.values[0]:12.8g}")


        # In[352]:


        #3. compute dc
        df_dynparam["dc_try"] = df_dynparam.apply(lambda x: (dc_min/dcmin_slip**p_dcscaleexp) * (x.slip_try)**p_dcscaleexp, axis=1)
        df_dynparam["hat_ds_factor_rapidnuc_nuc"] = hat_ds_factor_rapidnuc_nuc
        df_dynparam["hat_ds_factor_rapidnuc_patch"] = hat_ds_factor_rapidnuc_patch


        # In[353]:


        df_dynparam["dc_try"].sort_values() * 1e6


        # In[ ]:





        # In[354]:


        # plot debug Dc distribution
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        ax.plot(df_dynparam["slip_try"].values*1e6, df_dynparam["dc_try"].values*1e6, "o", c="k", alpha=0.5)

        # ax.text(0.02, 0.5, f"log slope={popt[1]:.3g}")

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel("Slip [μm]")
        ax.set_ylabel("Dc [μm]")

        # ax.set_xlim([1e-3, 1])
        # ax.set_ylim([1e-3, 1e-1])


        # In[355]:


        a_test = 2.4e-3
        a_test**3/a_patch**3


        # In[ ]:





        # In[356]:


        # plot debug for the estimation of minimum Dc and GIIC

        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        ax.plot(slip_all*1e6, GIIC_all, "o", c="gray", alpha=0.5)
        ax.plot(slip_syn*1e6, GIIC_syn, "k-")

        for index, df in df_dynparam.iterrows():
            ax.plot(np.array(df.slip_try)*1e6, np.array(0.5*df.dc_try*((df.fp_patch-df.hat_fr)*df.hat_sn_patch)), "o-", mfc="blue", mec="k", ms=5)
            
        ax.plot(dcmin_slip*1e6, np.array(0.5*dc_min*((df_dcmin.fp_patch-df_dcmin.hat_fr)*df_dcmin.hat_sn_patch)), "*", mfc="yellow", mec="k", ms=20)

        ax.text(0.02, 0.5, f"log slope={popt[1]:.3g}")
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel("Slip [μm]")
        ax.set_ylabel("Fracture energy [J/m$^2$]")

        ax.set_xlim([1e-3, 1])
        ax.set_ylim([1e-5, 0.5e1])

        fig.tight_layout()
        fig.savefig("../figure/debug_minDc_fitting.png", dpi=80)


        # ## Commpute the other model parameters

        # In[357]:


        def compute_GIIC(dc, fp, fr, sn):
            return 0.5*dc*(fp-fr)*sn


        # In[358]:


        # Peak and residual frictional resistance for the cohesive law
        # the friction coefficient is same between nuc and patch even for the rapid nucleation model
        df_dynparam.loc[:, "tau_c_nuc"] = df_dynparam.apply(lambda x: x.fp_patch*(hat_sn_patch*nuc_normalstress_alpha), axis=1) # same fp as the patch
        df_dynparam.loc[:, "tau_r_nuc"] = hat_fr_patch * (hat_sn_patch * nuc_normalstress_alpha) # same fr as the patch

        df_dynparam.loc[:, "tau_c_patch"] = df_dynparam.apply(lambda x: x.fp_patch*hat_sn_patch, axis=1)
        df_dynparam.loc[:, "tau_r_patch"] = hat_fr_patch*hat_sn_patch

        df_dynparam.loc[:, "tau_c_background"] = hat_fp_background * hat_sn_background
        df_dynparam.loc[:, "tau_r_background"] = hat_fr_background * hat_sn_background
                        
        # Initial stress state
        # df_dynparam.loc[:, "tau_0_nuc"] = df_dynparam.apply(lambda x: x.tau_c_patch*initialstress_fraction*(1+c_nucexcess), axis=1) # we control the stress in the nucleation area by gaussian distribution with 'c_nucexcess'
        df_dynparam.loc[:, "tau_0_nuc"] = df_dynparam.apply(lambda x: x.tau_c_patch*x.initialstress_fraction*(1+c_nucexcess), axis=1) # 
        # df_dynparam.loc[:, "tau_0_patch"] = df_dynparam.apply(lambda x: x.tau_c_patch*initialstress_fraction, axis=1)
        df_dynparam.loc[:, "tau_0_patch"] = df_dynparam.apply(lambda x: x.tau_c_patch*x.initialstress_fraction, axis=1)
        df_dynparam.loc[:, "tau_0_background"] = stressbackground_beta * hat_fr_background * hat_sn_background # for stress-free model, we investigate this value to arrest the rupture by the positive stress drop 

        # Compute fracture energy
        df_dynparam.loc[:, "GIIC_nuc"] = df_dynparam.apply(lambda x: compute_GIIC(x.dc_try, x.fp_patch, hat_fr_patch, (hat_sn_patch*nuc_normalstress_alpha)), axis=1) # Dc is same as the value in patch 
        df_dynparam.loc[:, "GIIC_patch"] = df_dynparam.apply(lambda x: compute_GIIC(x.dc_try, x.fp_patch, hat_fr_patch, hat_sn_patch), axis=1)
        df_dynparam.loc[:, "GIIC_background"] = 0.0



        # ## Compute critical nucleation size in nucleation area and patch region
        # We compute the critical nucleation size of 3D fault with the linear slip weakening law From Galis et al. (2014) as follows:
        # 
        # $$A_{init} = 1.75 S^{2.81} + 3.82 $$
        # 
        # $$ A_{init} = A_i / L_{fric}^2 $$
        # 
        # $$L_{fric} = \mu D_c / (\tau_s - \tau_d) $$

        # In[359]:


        def compute_Sratio(taup, taur, tau0):
            return (taup-tau0)/(tau0-taur)
            
        def compute_Ainit(S):
            return 1.75*(S**2.81) + 3.82

        def compute_Anuc(Ainit, mu, Dc, tptd):
            return Ainit*((mu*Dc) / (tptd))**2


        # In[360]:


        # compute Sratio
        df_dynparam.loc[:, "Sratio_nuc"] = 0 # assume tau0 = taup  as it is in critical state
        df_dynparam.loc[:, "Sratio_patch"] = df_dynparam.apply(lambda x: compute_Sratio(x.tau_c_patch, x.tau_r_patch, x.tau_0_patch), axis=1)

        # compute Ainit
        df_dynparam.loc[:, "Ainit_nuc"] =  df_dynparam.apply(lambda x: compute_Ainit(x.Sratio_nuc), axis=1)
        df_dynparam.loc[:, "Ainit_patch"] = df_dynparam.apply(lambda x: compute_Ainit(x.Sratio_patch), axis=1)

        # compute Anuc
        df_dynparam.loc[:, "Anuc_nuc"] =  df_dynparam.apply(lambda x: compute_Anuc(x.Ainit_nuc, mu, x.dc_try, (x.tau_c_nuc - x.tau_r_nuc)), axis=1)
        df_dynparam.loc[:, "Anuc_patch"] = df_dynparam.apply(lambda x: compute_Anuc(x.Ainit_patch, mu, x.dc_try, (x.tau_c_patch - x.tau_r_patch)), axis=1)

        # compute rnuc
        df_dynparam.loc[:, "rnuc_nuc"] =  df_dynparam.apply(lambda x: np.sqrt(x.Anuc_nuc/np.pi), axis=1) * 1e3
        df_dynparam.loc[:, "rnuc_patch"] = df_dynparam.apply(lambda x: np.sqrt(x.Anuc_patch/np.pi), axis=1) * 1e3


        # In[361]:


        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        df_dynparam.plot.scatter(x=["M0_mean"], y=["rnuc_patch"], ylabel="mm", ax=ax, label="patch")
        df_dynparam.plot.scatter(x=["M0_mean"], y=["rnuc_nuc"], ylabel="mm", ax=ax, c="r", label="nuc")
        ax.axhline(R_nuc*1e3, c="k", ls="--")

        ax.set_ylim([0, 4])


        # In[362]:


        gridsize_hpc = 0.041015625 #[mm]
        NgridperLc = df_dynparam.loc[129, "rnuc_patch"]/gridsize_hpc
        NgridperLc


        # In[ ]:





        # ## Check the amplitude of $f_p$

        # In[363]:


        df_dynparam_sorted = df_dynparam.sort_values("fp_patch").copy()

        fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))
        ax.plot(np.array(df_dynparam_sorted.M0_mean), np.array(df_dynparam_sorted.fp_patch), marker="o", c="k", ls="", ms=9, )

        # locs, _ = plt.xticks()
        # xlabels_all = [int(x.split('__')[1]) for x in df_modelparam.index]
        # plt.xticks(locs, xlabels_all)
        plt.axhline(1.0, c="k", ls="--")

        M0str = r"M$_{\mathrm{0}}$"

        ax.set_xlabel('{} [Nm]'.format(M0str))
        ax.set_ylabel("Peak friction coefficient $f_p$")

        ax.tick_params(axis='x', which='major', pad=5)

        ax.set_xlim([0, 1.4]);
        ax.set_ylim([0, 1.5]);

        ax.grid(True, which="major", c=np.array([230, 230, 230])/255, lw=0.25, zorder=-1)
        ax.set_axisbelow(True)
        fig.tight_layout()

        # plt.savefig(figdir + f"/peakfrictioncoef_allevents_{casename}.png", dpi=300, bbox_inches="tight")

        # ax.set_title("all the events")


        # In[ ]:





        # # Computational parameters

        # In[364]:


        # Some other model parameters related to the computation
        alpha_domain = 10 # 10 # computational domain size is alpha_domain times larger than radius of gauge
        nb_elements = 1024 #128 # #256 # number of grids per length; we set same for both x and z directions.

        # location of the center of nucleation zone
        nuc_z = 0

        # duration and time stepping
        if IfInitparam:
            duration = 2e-8 #8e-6 #10e-6 
        else:
            duration = 4.2e-6 #8e-6 #10e-6 dev for short duration
            
        tsf = 0.3 # factor of the critical time step

        char_reg_time = 0.0 #regularization time of normal stress for the case of linear coulomb friction law

        # dumping
        # dump_fields = "cohesion_0,cohesion_1,cohesion_2,top_disp_0,top_disp_1,top_disp_2,top_velo_0,top_velo_1,top_velo_2,G_c,tau_c,tau_r,mu_s,mu_k,d_c,load_0,load_1,load_2" # no space
        if IfInitparam:
            dump_fields = "cohesion_0,cohesion_1,top_disp_0,top_velo_0" # for HPC
        else:
            dump_fields = "cohesion_0,top_disp_0,top_velo_0" # for HPC

        # compute the total time step
        cs = np.sqrt(mu/rho)
        dx = float(R_patch*alpha_domain/nb_elements)
        dt_cfl = tsf*dx/cs
        Ntimestep = int(duration/dt_cfl)
        print(dt_cfl, Ntimestep)

        # we control the output frequency by the dt_dump
        dt_dump = 5e-8 # 2e-8 #1e-8 #2e-8 # [s] output sampling rate: dump data every dt_dump

        nb_dumps = int(np.ceil(Ntimestep/(dt_dump/dt_cfl)))
        nb_dumps, int(Ntimestep/nb_dumps)


        # In[365]:


        df_dynparam["fp_patch"]


        # In[366]:


        df_dynparam.columns


        # # Dump the file

        # In[ ]:





        # In[367]:


        # In[368]:


        # make list of parameters
        def generate_paramin_rapidnuc(df_model):
            param_in = []
            param_in.append(f"# The gouge patch dynamic rupture input file generated at {datetime.now().strftime('%Y-%m-%dT%H:%M:%S')}\n")
            param_in.append("string simulation_id = fb03-087__{:04d}_{:s}_{:.4f}\n".format(df_model.name, casename, delsigma_factor_range_event_case))
            
            param_in.append("\n# Computational domain size\n")
            param_in.append("double x_length = {:12.8e}\n".format(R_patch*alpha_domain))
            param_in.append("double z_length = {:12.8e}\n".format(R_patch*alpha_domain))
            
            param_in.append("int nb_x_elements = {:d}\n".format(nb_elements))
            param_in.append("int nb_z_elements = {:d}\n".format(nb_elements))
            
            # External loading
            param_in.append("\n# External loading\n")
            param_in.append("double sn_nuc = {:12.8e}\n".format(-hat_sn_patch*nuc_normalstress_alpha)) # sign convention is positive in opening
            param_in.append("double sn_patch = {:12.8e}\n".format(-hat_sn_patch))
            param_in.append("double sn_background = {:12.8e}\n".format(-hat_sn_background))
            
        #     param_in.append("double tau_nuc = {:12.8e}\n".format(df_model.tau_0_nuc))
            param_in.append("double tau_patch = {:12.8e}\n".format(df_model.tau_0_patch))
            param_in.append("double tau_background = {:12.8e}\n".format(df_model.tau_0_background))
            
            param_in.append("\n# Material constants\n")
            param_in.append("double E_top = {:12.8e}\n".format(E))
            param_in.append("double nu_top = {:12.8e}\n".format(nu))
            param_in.append("double rho_top = {:12.8e}\n".format(rho))
            
            # Frictional parameters; uniform GIIC in the patch
            param_in.append("\n# Frictional parameters\n")
            param_in.append("double Gc_nuc = {:12.8e}\n".format(df_model.GIIC_nuc))
            param_in.append("double tau_c_nuc = {:12.8e}\n".format(df_model.tau_c_nuc))
            param_in.append("double tau_r_nuc = {:12.8e}\n".format(df_model.tau_r_nuc))
            
            param_in.append("double Gc_patch = {:12.8e}\n".format(df_model.GIIC_patch))
            param_in.append("double tau_c_patch = {:12.8e}\n".format(df_model.tau_c_patch))
            param_in.append("double tau_r_patch = {:12.8e}\n".format(df_model.tau_r_patch))
            
        #     param_in.append("double Gc_background = {:12.8e}\n".format(df_model.GIIC_patch)) # To avoid the zero division, we set the same GIIC as the patch; not used in the simulation with Coulomb's law
            param_in.append("double tau_c_background = {:12.8e}\n".format(df_model.tau_c_background))
            param_in.append("double tau_r_background = {:12.8e}\n".format(df_model.tau_r_background))
            
            param_in.append("double dc_nuc = {:12.8e}\n".format(df_model.dc_try)) # same as the patch
            param_in.append("double ds_nuc = {:12.8e}\n".format(hat_ds_factor_rapidnuc_nuc * df_model.dc_try))
        #     param_in.append("double ds_nuc = {:12.8e}\n".format(hat_ds_nuc))
            param_in.append("double fp_nuc = {:12.8e}\n".format(df_model.fp_patch)) # same as the fp patch, which is modified when using fp Gaussian nucleation
        #     param_in.append("double fr_nuc = {:12.8e}\n".format(df_model.hat_fr/nuc_normalstress_alpha)) # after the nucleation, stress drop is same level with the patch domain #(df_model.hat_fr))
            param_in.append("double fr_nuc = {:12.8e}\n".format(df_model.hat_fr)) # set same fr as the patch
        
            param_in.append("double dc_patch = {:12.8e}\n".format(df_model.dc_try))
            param_in.append("double ds_patch = {:12.8e}\n".format(hat_ds_factor_rapidnuc_patch * df_model.dc_try))
        #     param_in.append("double ds_patch = {:12.8e}\n".format(hat_ds_patch))
            param_in.append("double fp_patch = {:12.8e}\n".format(df_model.fp_patch))
            param_in.append("double fr_patch = {:12.8e}\n".format(df_model.hat_fr))
            
            param_in.append("double dc_background = {:12.8e}\n".format(df_model.dc_try)) # need to be non-zero to avoid zero devision
            param_in.append("double ds_background = {:12.8e}\n".format(hat_ds_factor_rapidnuc_background * df_model.dc_try))
            # param_in.append("double dc_background = {:12.8e}\n".format(100*df_model.dc_try)) # need to be non-zero to avoid zero devision
            # param_in.append("double ds_background = {:12.8e}\n".format(100*hat_ds_factor_rapidnuc_background * df_model.dc_try))
            param_in.append("double fp_background = {:12.8e}\n".format(df_model.hat_fp_background)) #### peak is same as redisual in background
            param_in.append("double fr_background = {:12.8e}\n".format(df_model.hat_fr_background))
            

            # patch sizes
            param_in.append("\n# Patch sizes and locations\n")
            param_in.append("double R_nuc = {:12.8e}\n".format(R_nuc)) # sign convention is positive in opening
            param_in.append("double R_patch = {:12.8e}\n".format(R_patch))
            param_in.append("double R_margin = {:12.8e}\n".format(R_margin))
            
            param_in.append("double nuc_x = {:12.8e}\n".format(nuc_x))
            param_in.append("double nuc_z = {:12.8e}\n".format(nuc_z))

            param_in.append("double c_nucexcess = {:12.8e}\n".format(c_nucexcess))
        #     param_in.append("double alpha_fpnuc = {:12.8e}\n".format(alpha_fpnuc)) 
        #     param_in.append("double c_fpnuc = {:12.8e}\n".format(c_fpnuc))

            # simulation parameters
            param_in.append("\n# Simulation parameters\n")
            param_in.append("double duration = {:12.8e}\n".format(duration))
            param_in.append("double tsf = {:12.8e}\n".format(tsf))
            param_in.append("double char_reg_time= {:12.8e}\n".format(char_reg_time))
            
            param_in.append("string dump_fields = {:s}\n".format(dump_fields))
            param_in.append("int nb_dumps = {:d}\n".format(nb_dumps))

            return param_in


        # ## Write the input files

        # In[369]:


        # select some cases for preliminary result
        # selectids = [24, 50, 52, 72, 129] # New set of non-self-similar events with merged event catalog.
        # selectids = [50, 52, 72, 129] # New set of non-self-similar events with merged event catalog.
        # selectids = [50, 52] # New set of non-self-similar events with merged event catalog.
        # selectids = [24] # New set of non-self-similar events with merged event catalog.



        # In[370]:


        for i, df_model in df_dynparam.iterrows():
            # if i not in selectids:
            if i != selectid:
                # print(i)
                continue

            if IfInitparam:
                foname = inputfileoutdir+f"/rupgougepatch_fb03-{expr_id:03d}__{df_model.name:04d}_{casename}_{delsigma_factor_range_event_case:.4f}_initparam.in"
            else:
                foname = inputfileoutdir+f"/rupgougepatch_fb03-{expr_id:03d}__{df_model.name:04d}_{casename}_{delsigma_factor_range_event_case:.4f}.in"

            print(f"output {foname}")
        #     param_in = generate_paramin(df_model)
        #     param_in = generate_paramin_linear_coulomb_friction_law(df_model)
            param_in = generate_paramin_rapidnuc(df_model)
            # output file
            with open(foname, "w") as fo:
                fo.writelines(param_in)


        # In[371]:


        # Dump the csv file
        # df_dynparam.to_csv(f"../data/gouge_dynamicrupture_modelparam_{casename}.csv", float_format="%12.8e")


        # In[ ]:





        # In[ ]:




