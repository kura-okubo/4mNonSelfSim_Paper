#!/usr/bin/env python

# Functions used to conduct MCMC inversion of AE source properties.
# 2022.05.09 Kurama Okubo
# 2024.02.07 update for aperture effect correction
# 2024.06.09 update for Küpper wavelet

import os
import obspy
from obspy import read, Stream, Trace
from scipy import signal
import matplotlib.pyplot as plt
import glob
from glob import glob
import numpy as np
import mpmath as mp
import pandas as pd
import datetime
# from datetime import timedelta # as datetime.time delta is not working in nanosecond, use UTCDateTime
from tqdm import tqdm
import pickle
import warnings
from obspy.core.utcdatetime import UTCDateTime    
import time
from scipy.spatial.distance import euclidean, correlation
from obspy.signal.cross_correlation import correlate, xcorr_max
from scipy.optimize import fsolve

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 12
os.environ['TZ'] = 'GMT' # change time zone to avoid confusion in unix_tvec conversion
UTCDateTime.DEFAULT_PRECISION = 8 # increase the time precision


def M02Mw(M0):
    """
    convert from M0 to Mw
    """
    return (np.log10(M0) - 9.105) * 2.0 / 3.0

def Mw2M0(Mw):
    """
    convert from Mw to M0
    """
    return 10**( 1.5 * Mw + 9.105)

def stf_cosine(t, TR, fz):
    '''
        source time function of cosine wavelet
        https://openswpc.github.io/2._Parameters/0207_source/
        Argument: 
            t:: time vector
    '''
    stf = np.zeros(len(t))
    for i, tt in enumerate(t):
        if 0<tt and tt<TR:
            stf[i] = (fz/TR) * (1 - np.cos((2*np.pi*tt/TR)))
        else:
            stf[i] = 0
            
    return stf

def stf_kupper(t, TR, fz):
    '''
        source time function of Küpper wavelet
        https://openswpc.github.io/2._Parameters/0207_source/
        Argument: 
            t:: time vector
    '''
#     print("kupper test")
    stf = np.zeros(len(t))
    for i, tt in enumerate(t):
        if 0<tt and tt<TR:
            stf[i] = fz * ((3 * np.pi)/(4 * TR)) * np.sin(np.pi * tt/TR)**3
        else:
            stf[i] = 0
            
    return stf

def compute_momenttensor_lab(strike_deg, rake_deg):
    '''
    Return moment tensor associated with events on 4m biax fault
    assuming that dip = 90 (sidecoord: x along fault, y along side of fault, z upwards.) 
    input:
        rake (radian): rake of slip. Positive rotating from x to z direction.
    output:
        mij: moment tensor normalized such that np.linalg.norm(mij, ord=None) = 1 following OpenSWPC source normalization (Frobenius norm)
    '''
    
    strike = np.deg2rad(strike_deg)
    rake = np.deg2rad(rake_deg)
    
    mij = np.zeros((3, 3))
    mij[0][0] = -np.cos(rake) * np.sin(2*strike)
    mij[1][1] = np.cos(rake) * np.sin(2*strike)
    mij[2][2] = 0
    mij[0][1] = np.cos(rake) * np.cos(2*strike)
    mij[0][2] = -np.sin(rake) * np.sin(strike)
    mij[1][0] = mij[0][1]
    mij[1][2] = np.sin(rake)*np.cos(strike) # positive in the rotation from x to y direction. 
    mij[2][0] = mij[0][2]
    mij[2][1] = mij[1][2]
    
    mag = np.linalg.norm(mij, ord=None)
    
    return mij/mag

def compute_biax_synthesize_waveform_vz(M0hat, rake_deg, TR_new, TR_input, st_syn_sta, k=1e-6, strike_deg=0):
    """
    synthesize the waveform with 3 source parameters.
    Input:
    M0hat: seismic moment multiplied by sqrt(2); M0hat=sqrt(2)*M0 [Nm]
    rake_deg: rake of slip [deg]
    TR_new: new source duration [s] used to estimate dominant frequency.
    TR_input : source duration [s] precompiled in the OpenSWPC.
    st_syn_sta: Stream of Green's function convolved be pre-convolved STF.
    k: coefficient of waterlevel stabilization (default: 1.0e-6)
    strike_deg: Strike (incliment of fault) [deg] (default: 0 i.e. horizontal)
    
    Workflow:
    1. reconvolve STF on the Green's functions independently 
        as each Green's function is convolved with the input STF defined in the OpenSWPC input file.
    2. compute mij with rake assuming double couple source on the fault.
    3. sum the components to synthesize the waveform.
    
    Note: to syncronize the coordinate system on the isocoord used in the numerical simulation,
        flip the sign associated with x and z direction.
        The sensors on southside (>=OL17) is mapped across the fault. Then, compute the Green's function on isocoord,
        where x and z axis is flipped from the absolute coordinate system. The resulting Green's function is associated with
        isocomp/double couple source on the x-and-z-flipped coordinate system.
        Thus, we need to flip back to convert to the absolute coordinate system on the Mxx, Mxy, Mzz, Myz.
        Then, vz is defined as positive in the direction of normal to the side surface for both north and south side of sensors.
        The preamp is negative filter, while the compression of sensor surface is negative.
        Thus the observed signal is also positive in the direction of normal of side surfaces.
    """
    
    #----------------------------------#
    #1. reconvolve new STF---#
    #----------------------------------#
    
    dt = st_syn_sta[0].stats.delta
    
    # Source time function pre-convolved in the OpenSWPC
    tvec_STF_input = np.arange(np.ceil(TR_input/dt)) * dt # time vector of pre-convolved STF [J/s]
    # stf_cos_input = stf_cosine(tvec_STF_input, TR_input, 1.0) # the STF pre-convolved in the OpenSWPC is scaed to be unity.
    stf_kupper_input = stf_kupper(tvec_STF_input, TR_input, 1.0) # the STF pre-convolved in the OpenSWPC is scaed to be unity.
    
    # Source time function of new target STF
    tvec_STF_new = np.arange(np.ceil(TR_new/dt)) * dt # time vector of new STF
    # stf_cos_new = stf_cosine(tvec_STF_new, TR_new, 1.0) # reconvolve cosine source with TR_new
    stf_kupper_new = stf_kupper(tvec_STF_new, TR_new, 1.0) # reconvolve cosine source with TR_new

    # check if Green's tensor comprises of six components
    if len(st_syn_sta) != 6:
        warnings.warn("number of green's tensor is not 6.")
  
    st_syn_reconv = st_syn_sta.copy() # copy stream containing the reconvolved Green's function

    #---reconvolve STF---#
    N = st_syn_reconv[0].stats.npts # number of points on the green's function data
    FS0 = np.fft.rfft(stf_kupper_input, N) # FFT of pre-convolved STF
    FS1 = np.fft.rfft(stf_kupper_new, N) # FFT of new STF
    tmp1 = FS1 * np.conjugate(FS0) # dot product
    tmp2 = np.maximum(np.abs(FS0)**2, (k*np.abs(FS0).max())**2)
    for i, tr in enumerate(st_syn_reconv):
        FX0 = np.fft.rfft(tr.copy().taper(0.05).data, N)
        FX1= (tmp1/tmp2)*FX0
        st_syn_reconv[i].data = np.fft.irfft(FX1, N).real
    #--------------------------#
    
    Gzxx = st_syn_reconv.select(channel="G_Vz_mxx")[0].data
    Gzyy = st_syn_reconv.select(channel="G_Vz_myy")[0].data
    Gzzz = st_syn_reconv.select(channel="G_Vz_mzz")[0].data
    Gzxy = st_syn_reconv.select(channel="G_Vz_mxy")[0].data
    Gzxz = st_syn_reconv.select(channel="G_Vz_mxz")[0].data
    Gzyz = st_syn_reconv.select(channel="G_Vz_myz")[0].data
    
    #----------------------------------#
    #2. compute mij ------------#
    #----------------------------------#
    mij0 = compute_momenttensor_lab(strike_deg, rake_deg)
    
    #----------------------------------#
    #3. synthesize waveform -#
    #----------------------------------#
    
    # convert from isocoord to absolute experimental coordinate system
    # i.e. flip on x and z axis
    station_id = int(st_syn_sta[0].stats.station.split('OL')[-1])
    
    if station_id >= 17: # sensors installed on the north side of rock
        Gzxx = -Gzxx
        Gzzz = -Gzzz
        Gzxy = -Gzxy
        Gzyz = -Gzyz
    
    # sum up the moment tensors
    vz = M0hat * (mij0[0, 0] * Gzxx + mij0[1, 1] * Gzyy + mij0[2, 2] * Gzzz +
               mij0[0, 1] * Gzxy + mij0[0, 2] * Gzxz + mij0[1, 2] * Gzyz)

    return vz, FS0, FS1, tmp1, tmp2


#---Functions of Gridsearch---$

# Compute amplitude response with incident angle
def incidentangle_scalingfactor_analytic(v, theta, TR, R):
    if theta==0:
        return 1.0
    else:
        va = v/np.sin(theta)
        J1 = mp.besselj(1, (2*np.pi*R)/(va*TR))
        return  ((va * TR)/(np.pi*R)) * J1
    
def get_retrimmed_traces(M0hat_try, rake_try, TR_try, param, fix_dtshift=False, aperturecorrection=False):
    """
        synthesize the traces with given M0hat, rake_deg and TR.
        correct the arrival time using cross-correlation, and retrim the traces
        return the all traces. We detached the VR computation as different function. 
    """

    # initialize paramters
    st_event = param["st_event"].copy()
    TR_input = param["TR_input"] # set at the input file of OpenSWPC
    cp = param["cp"] #[m/s] dilational wave velocity
    cs = param["cs"] #[m/s] shear wave velocity
    VR_sensors = param["VR_sensors"] # list of sensor to compute VR

    prePwinlen = param["prePwinlen"] # [ms] window length behind the p wave arrival
    Pwinlen = param["Pwinlen"] # [ms] window length behind the p wave arrival
    preSwinlen = param["preSwinlen"] # [ms] window length behind the p wave arrival
    Swinlen = param["Swinlen"] # [ms] window length behind the p wave arrival

    ifFilteringBeforeTrim = param["ifFilteringBeforeTrim"]
    weight_factor = param["weight_factor"]
    freqmin = param["freqmin"] #[Hz]
    freqmax = param["freqmax"] #[Hz]
    max_lag_shift = param["max_lag_shift"] # number of time shift for cross-correlation

    UTCDateTime.DEFAULT_PRECISION = 9 # need to correctly compute in ns order


    sensorids = np.unique([tr.stats.station for tr in st_event])

    # start grid search

    #--------------------------------------------------------------------#
    #---Synthesize the waveform and store to the stream--#
    #--------------------------------------------------------------------#
    # adjust k_waterlevel with respect to TR_new
    if TR_try <= 5e-6:
        k_watarlevel = 1e-3 # waterlevel
    else:
        k_watarlevel = 1e-4 # waterlevel

    st_retrim = Stream()
    
    for sensorid in VR_sensors: # compute the component of VR within VR_sensors
        st_syn_sta = st_event.select(station=sensorid, location='GF')
        # applying taper on the synthetic green's functions
        st_syn_sta.taper(0.05)

        # st_syn_sta.detrend(type='linear')
        # st_syn_sta.filter("lowpass", freq = 0.6e6, corners=2, zerophase=True)

        vz, _, _, _, _ = compute_biax_synthesize_waveform_vz(M0hat_try, rake_try, TR_try, TR_input, st_syn_sta, k=k_watarlevel, strike_deg=0)
        
        # store to the event Stream
        tr_vz = st_syn_sta[0].copy()
        tr_vz.data = vz
        tr_vz.stats.location = "Syn"
        tr_vz.stats.channel = "VZ"

        st_event.append(tr_vz)

        # correct arrival time with cross-correlation and compute trimmed trace.
#         for sensorid in sensorids:
        tr_obs = st_event.select(station=sensorid, location="stage1", channel="OZ")[0] # select the response removed traces
        tr_syn = st_event.select(station=sensorid, channel="VZ")[0]
        p_arrival = tr_obs.stats.dist/cp #[ms]
        s_arrival = tr_obs.stats.dist/cs #[ms]


        # Apply filter before trimming and computing the residual
        if ifFilteringBeforeTrim:
            tr_obs_filtered = tr_obs.copy().taper(0.05).filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
            tr_syn_filtered = tr_syn.copy().taper(0.05).filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
        else:
            tr_obs_filtered = tr_obs.copy().taper(0.05)
            tr_syn_filtered = tr_syn.copy().taper(0.05)


        tr_obs_filtered.stats.location = "Filt"
        tr_syn_filtered.stats.location = "Filt"

        # store p and s arrival time with pretrigger
        tr_obs_filtered.stats.p_arrival = p_arrival
        tr_obs_filtered.stats.s_arrival = s_arrival

        st_event.append(tr_obs_filtered)
        st_event.append(tr_syn_filtered)

        #--------------------------------------------------------------------#
        #---Correction of arrival time by cross-correlation--------#
        # We first temporally trim the window, then compute cc to find the shift to maximize the cc.
        # We then correct the start time on synthetic trace, then re-trimming the trance.
        #--------------------------------------------------------------------#   

        #------------------------------#
        #---process on P wave window---#
        #------------------------------#

        st_tmp = tr_syn.stats.starttime
        # pt = timedelta(milliseconds = tr_obs.stats.pretrigger)
        # pp = timedelta(milliseconds = p_arrival)
        pt = tr_obs.stats.pretrigger * 1e-3 #[s]
        pp = p_arrival * 1e-3 #[s]
        # residu_init_p =  timedelta(milliseconds = p_arrival-prePwinlen)
        residu_init_p =  (p_arrival-prePwinlen) * 1e-3 #[s]

        # residu_winlen_only_P = timedelta(milliseconds = Pwinlen) # P window with fixed length
        residu_winlen_only_P =Pwinlen * 1e-3 # [s] P window with fixed length

        starttime_p = st_tmp+(pt+residu_init_p) # plus with [s] to UTCtimedate
        endtime_p  = st_tmp+(pt+pp+residu_winlen_only_P) # plus with [s] to UTCtimedate
        # print("Pwinlen", pt, pp, residu_init_p, residu_winlen_only_P, starttime_p, endtime_p)
#             print(starttime_tmp, endtime_tmp)

        # Trim the trial traces
        tr_obs_trim_p = tr_obs_filtered.copy().trim(starttime_p, endtime_p, pad=True, fill_value=0, nearest_sample=True)
        tr_syn_trim_p = tr_syn_filtered.copy().trim(starttime_p, endtime_p, pad=True, fill_value=0, nearest_sample=True)

        if fix_dtshift:
            # print("debug")
            dt_shift = param["fix_dtshift"][sensorid][0]
            cc=0
            N_shift=int(dt_shift/tr_obs_trim_p.stats.delta)

        else:
            # Compute cross-correlation
            cc = correlate(tr_obs_trim_p.taper(0.1), tr_syn_trim_p.taper(0.1), max_lag_shift, demean=True, normalize='naive') # the order is 1. obs and 2. syn
            N_shift, _ = xcorr_max(cc,  abs_max=False)
            # print(cc[0], N_shift)
            dt_shift = N_shift * tr_obs_trim_p.stats.delta # [s] time to shift for the correction of arrival
        
        # Shift the start time of synthetic filtered waveform
        tr_syn_filtered_shifted_p = tr_syn_filtered.copy()
        tr_syn_filtered_shifted_p.stats.starttime += dt_shift # shift the trace; the order of nanosecond is correctly stored into starttime.
        tr_syn_filtered_shifted_p.stats.location = "FiltShifted_P"
        tr_syn_filtered_shifted_p.stats.cc = cc
        tr_syn_filtered_shifted_p.stats.N_shift = N_shift
        tr_syn_filtered_shifted_p.stats.dt_shift = dt_shift
        st_event.append(tr_syn_filtered_shifted_p)

        # retrim p window with shifted p winlen
        # pp_shift = timedelta(milliseconds = dt_shift*1e3)
        pp_shift = dt_shift #[s]
        starttime_p_shift = st_tmp+(pt+residu_init_p + pp_shift) # plus with [s] to UTCtimedate
        endtime_p_shift  = st_tmp+(pt+pp+residu_winlen_only_P + pp_shift) # plus with [s] to UTCtimedate

        tr_obs_retrimmed_p = tr_obs_filtered.copy().trim(starttime_p_shift, endtime_p_shift, pad=True, fill_value=0, nearest_sample=True)            
        tr_syn_retrimmed_p = tr_syn_filtered_shifted_p.copy().trim(starttime_p_shift, endtime_p_shift, pad=True, fill_value=0, nearest_sample=True)

        if aperturecorrection:
            # divided the amplitude by the factor of aperture effect
            df_incidentangle_sensor = param["df_incidentangle"].loc[f"{sensorid}__{param['datacase']}"]
            beta_coef_p = float(incidentangle_scalingfactor_analytic(param["cp"], np.deg2rad(df_incidentangle_sensor["incidentangle"]), TR_try, param["R_sensor"]))
            print(f"aperture correct: {sensorid} incident angle={df_incidentangle_sensor['incidentangle']:4f}deg. beta_p={beta_coef_p:4f}")
            tr_obs_retrimmed_p.data /= beta_coef_p # correct the observed P wave amplitude
            
        tr_obs_retrimmed_p.stats.location = "Retrim_P"
        tr_obs_retrimmed_p.stats.starttime_p = starttime_p
        tr_obs_retrimmed_p.stats.endtime_p = endtime_p
        tr_obs_retrimmed_p.stats.starttime_p_shift = starttime_p_shift
        tr_obs_retrimmed_p.stats.endtime_p_shift = endtime_p_shift

        tr_obs_retrimmed_p.stats.pwin_start_tmp = residu_init_p # p tmp window starttime from pretrigger
        tr_obs_retrimmed_p.stats.pwin_end_tmp = pp+residu_winlen_only_P # p tmp window endtime from pretrigger

        tr_obs_retrimmed_p.stats.pwin_start = residu_init_p + pp_shift # p window starttime from pretrigger
        tr_obs_retrimmed_p.stats.pwin_end = pp+residu_winlen_only_P +pp_shift # p window endtime from pretrigger

        tr_syn_retrimmed_p.stats.location = "Retrim_P"
        tr_syn_retrimmed_p.stats.cc = cc
        tr_syn_retrimmed_p.stats.N_shift = N_shift
        tr_syn_retrimmed_p.stats.dt_shift = dt_shift
        
        # print(N_shift)
        # print(pp_shift)
        # print(tr_obs_filtered)
        # print(tr_syn_filtered_shifted_p)
        # print(tr_obs_retrimmed_p)
        # print(tr_syn_retrimmed_p)

        # Update: set the number of the data sample as the same in the case when the trim window causes one data point difference
        if tr_obs_retrimmed_p.stats.npts > tr_syn_retrimmed_p.stats.npts:
            Ndiff = tr_obs_retrimmed_p.stats.npts-tr_syn_retrimmed_p.stats.npts
            tr_obs_retrimmed_p.data = tr_obs_retrimmed_p.data[:-Ndiff]
        elif tr_obs_retrimmed_p.stats.npts < tr_syn_retrimmed_p.stats.npts:
            Ndiff = tr_syn_retrimmed_p.stats.npts-tr_obs_retrimmed_p.stats.npts
            tr_syn_retrimmed_p.data = tr_syn_retrimmed_p.data[:-Ndiff]

        assert tr_obs_retrimmed_p.stats.npts == tr_syn_retrimmed_p.stats.npts # check if traces have the same data point

        tr_obs_trim_p.stats.location = "TrimTmp_P"
        tr_syn_trim_p.stats.location = "TrimTmp_P"
        # save the trimming time

        st_event.append(tr_obs_trim_p)
        st_event.append(tr_syn_trim_p)
        st_event.append(tr_obs_retrimmed_p)
        st_event.append(tr_syn_retrimmed_p)

        #------------------------------#
        #---process on S wave window---#
        #------------------------------#

        # pt = timedelta(milliseconds = tr_obs.stats.pretrigger)
        # ps = timedelta(milliseconds = s_arrival)
        # residu_init_s =  timedelta(milliseconds = s_arrival-preSwinlen)

        pt = tr_obs.stats.pretrigger * 1e-3 #[s]
        ps = s_arrival * 1e-3 #[s]
        residu_init_s =  (s_arrival-preSwinlen) * 1e-3 #[s]

        residu_winlen_only_S = Swinlen * 1e-3 # [s] P window with fixed length

        starttime_s = st_tmp+(pt+residu_init_s) # plus with [s] to UTCtimedate
        endtime_s  = st_tmp+(pt+ps+residu_winlen_only_S) # plus with [s] to UTCtimedate
        # print("Swinlen", pt, ps, residu_init_s, residu_winlen_only_S)
#             print(starttime_tmp, endtime_tmp)

        # Trim the trial traces
        tr_obs_trim_s = tr_obs_filtered.copy().trim(starttime_s, endtime_s, pad=True, fill_value=0, nearest_sample=True)
        tr_syn_trim_s = tr_syn_filtered.copy().trim(starttime_s, endtime_s, pad=True, fill_value=0, nearest_sample=True)

        # Compute cross-correlation
        if fix_dtshift:
            dt_shift = param["fix_dtshift"][sensorid][1]
            cc=0
            N_shift=int(dt_shift/tr_obs_trim_s.stats.delta)

        else:
            cc = correlate(tr_obs_trim_s.taper(0.1), tr_syn_trim_s.taper(0.1), max_lag_shift, demean=True, normalize='naive') # the order is 1. obs and 2. syn
            N_shift, _ = xcorr_max(cc,  abs_max=False)
            # print(cc[0], N_shift)
            dt_shift = N_shift * tr_obs_trim_s.stats.delta # [s] time to shift for the correction of arrival
        
        # Shift the start time of synthetic filtered waveform
        tr_syn_filtered_shifted_s = tr_syn_filtered.copy()
        tr_syn_filtered_shifted_s.stats.starttime += dt_shift # shift the trace; the order of nanosecond is correctly stored into starttime.
        tr_syn_filtered_shifted_s.stats.location = "FiltShifted_S"
        tr_syn_filtered_shifted_s.stats.cc = cc
        tr_syn_filtered_shifted_s.stats.N_shift = N_shift
        tr_syn_filtered_shifted_s.stats.dt_shift = dt_shift
        st_event.append(tr_syn_filtered_shifted_s)

        # retrim s window with shifted s winlen
        # ps_shift = timedelta(milliseconds = dt_shift*1e3)
        ps_shift = dt_shift
        starttime_s_shift = st_tmp+(pt+residu_init_s + ps_shift) # plus with [s] to UTCtimedate
        endtime_s_shift  = st_tmp+(pt+ps+residu_winlen_only_S + ps_shift) # plus with [s] to UTCtimedate
        tr_obs_retrimmed_s = tr_obs_filtered.copy().trim(starttime_s_shift, endtime_s_shift, pad=True, fill_value=0, nearest_sample=True)            
        tr_syn_retrimmed_s = tr_syn_filtered_shifted_s.copy().trim(starttime_s_shift, endtime_s_shift, pad=True, fill_value=0, nearest_sample=True)

        if aperturecorrection:
            # divided the amplitude by the factor of aperture effect
            beta_coef_s = float(incidentangle_scalingfactor_analytic(param["cs"], np.deg2rad(df_incidentangle_sensor["incidentangle"]), TR_try, param["R_sensor"]))
            # print(f"aperture correct: {sensorid} incident angle={df_incidentangle_sensor['incidentangle']:4f}deg. beta_s={beta_coef_s:4f}")
#             tr_obs_retrimmed_s.data /= beta_coef_s # correct the observed P wave amplitude: The S correction does not work well. Skipping.
            
        tr_obs_retrimmed_s.stats.location = "Retrim_S"
        tr_obs_retrimmed_s.stats.starttime_s = starttime_s
        tr_obs_retrimmed_s.stats.endtime_s = endtime_s
        tr_obs_retrimmed_s.stats.starttime_s_shift = starttime_s_shift
        tr_obs_retrimmed_s.stats.endtime_s_shift = endtime_s_shift
        tr_obs_retrimmed_s.stats.swin_start_tmp = residu_init_s# p window starttime from pretrigger
        tr_obs_retrimmed_s.stats.swin_end_tmp = ps+residu_winlen_only_S # p window endtime from pretrigger
        tr_obs_retrimmed_s.stats.swin_start = residu_init_s + ps_shift# p window starttime from pretrigger
        tr_obs_retrimmed_s.stats.swin_end = ps+residu_winlen_only_S + ps_shift # p window endtime from pretrigger

        tr_syn_retrimmed_s.stats.location = "Retrim_S"
        tr_syn_retrimmed_s.stats.cc = cc
        tr_syn_retrimmed_s.stats.N_shift = N_shift
        tr_syn_retrimmed_s.stats.dt_shift = dt_shift
        
        # print(N_shift)
        # print(ps_shift)
        # print(tr_obs_filtered)
        # print(tr_syn_filtered_shifted_s)
        # print(tr_obs_retrimmed_s)
        # print(tr_syn_retrimmed_s)
        # print(dt_shift)

        # set the number of the data sample as the same
        if tr_obs_retrimmed_s.stats.npts > tr_syn_retrimmed_s.stats.npts:
            Ndiff = tr_obs_retrimmed_s.stats.npts-tr_syn_retrimmed_s.stats.npts
            tr_obs_retrimmed_s.data = tr_obs_retrimmed_s.data[:-Ndiff]
        elif tr_obs_retrimmed_s.stats.npts < tr_syn_retrimmed_s.stats.npts:
            Ndiff = tr_syn_retrimmed_s.stats.npts-tr_obs_retrimmed_s.stats.npts
            tr_syn_retrimmed_s.data = tr_syn_retrimmed_s.data[:-Ndiff]

        # print( tr_obs_retrimmed_s.stats.npts, tr_syn_retrimmed_s.stats.npts )
        assert tr_obs_retrimmed_s.stats.npts == tr_syn_retrimmed_s.stats.npts # check if traces have the same data point

        tr_obs_trim_s.stats.location = "TrimTmp_S"
        tr_syn_trim_s.stats.location = "TrimTmp_S"
        # save the trimming time

        st_event.append(tr_obs_trim_s)
        st_event.append(tr_syn_trim_s)
        st_event.append(tr_obs_retrimmed_s)
        st_event.append(tr_syn_retrimmed_s)


    # #--------------------------------------------------------------------#
    # #---Compute the variance reduction--#
    # #--------------------------------------------------------------------#   
    # VR_sensor = {} # dictionary to store VR of each sensor

    # for sensorid in sensorids:
        
    #     tr_obs_retrimmed = st_event.select(station=sensorid, location="Retrim", channel="OZ")[0]
    #     tr_syn_retrimmed = st_event.select(station=sensorid, location="Retrim", channel="VZ")[0]
        
    #     VRtmp1 = np.linalg.norm(tr_obs_retrimmed.data - 1.0 * tr_syn_retrimmed.data, axis=0) ** 2 # here the synthetic has been already scaled
    #     VRtmp2 = np.linalg.norm(tr_obs_retrimmed.data, axis=0)**2
    #     VR_sensor[sensorid] = weight_all[sensorid] * (VRtmp1/VRtmp2)

    # # sum up to compute VR
    # VR = 1 - sum(VR_sensor.values())
    
    return st_event #, weight_all
 
def compute_p_reflection(tr_obs, param):
    cp = param["cp"] #[m/s]
    rock_width = 0.1 #[m] 100mm width of fault
    rock_hight = 0.2 #[m] 200mm on onside of rock specimen
    sensorloc = np.array([0, 0, 70e-3])
    sourceloc = np.array([tr_obs.stats.xi1, tr_obs.stats.eta1, tr_obs.stats.zeta1])

    rvec = sourceloc - sensorloc
    wvec_side = [0, -2*rock_width, 0]
    wvec_bottom = [0, 0, 2*(rock_hight-sensorloc[2])]
    wvec_top = [0, 0, -2*(rock_hight+sensorloc[2])]

    def f_lawcos(x):
        return np.sqrt( np.linalg.norm(rvec)**2 + np.linalg.norm(x)**2 - 2 * np.dot(rvec, x))

    # direct p
    p_direct = np.linalg.norm(rvec) / cp
    p_direct_ref =  f_lawcos([0, 0, 0]) / cp
    assert p_direct == p_direct_ref

    # pp from side
    pp_side = f_lawcos(wvec_side) / cp

    # pp from bottom
    pp_bottom = f_lawcos(wvec_bottom) / cp

    # pp from top
    pp_top = f_lawcos(wvec_top) / cp

    #ppp from side
    ppp_side = np.sqrt( np.linalg.norm(rvec)**2 + np.linalg.norm(-np.array(wvec_side))**2 - 2 * np.dot(rvec, -np.array(wvec_side))) / cp
    # print(pp_side, ppp_side)

    return (p_direct, pp_side, pp_bottom, pp_top, ppp_side)



def compute_ps_and_sp_side_reflection(tr_obs, param):
    """
    compute ps and sp reflection from side surface by solving the simultaneous equation with Snell's law 
    """
    cp = param["cp"] #[m/s]
    cs = param["cs"] #[m/s]

    rock_width = 0.1 #[m] 100mm width of fault
    rock_hight = 0.2 #[m] 200mm on onside of rock specimen
    sensorloc = np.array([0, 0, 70e-3])
    sourceloc = np.array([tr_obs.stats.xi1, tr_obs.stats.eta1, tr_obs.stats.zeta1])

    rvec = np.array(sourceloc - sensorloc)
    wvec_side = np.array([0, -2*rock_width, 0])

    # Function to compute PS
    # define m vector from sensor to reflection point
    vnorm = np.linalg.norm

    def func_ps(x):
        p, q = x
        if p<0 or q<0:
            return [np.inf, np.inf] 
        
        mvec = p*rvec + q*wvec_side
        svec = mvec - rvec
        cos_ps_i = np.dot(svec, wvec_side)/(vnorm(svec) * vnorm(wvec_side))
        cos_ps_j = np.dot(mvec, wvec_side)/(vnorm(mvec) * vnorm(wvec_side))
        sin_ps_i = np.sqrt(1-cos_ps_i**2)
        sin_ps_j = np.sqrt(1-cos_ps_j**2)
        f = np.zeros(2)
        f[0] = np.dot(mvec, wvec_side) / vnorm(wvec_side)**2 - 0.5
        f[1] = sin_ps_i/cp - sin_ps_j/cs
        return f

    A = fsolve(func_ps, [0.1, 0.1], full_output=True, xtol=1.49012e-08)
    print(A)
    # compute i and j
    p0_ps, q0_ps = A[0]
    mvec_ps_side = p0_ps*rvec + q0_ps*wvec_side
    svec_ps_side = mvec_ps_side - rvec

    ideg_ps = np.rad2deg(np.arccos(np.dot(svec_ps_side, wvec_side)/(vnorm(svec_ps_side) * vnorm(wvec_side))))
    jdeg_ps = np.rad2deg(np.arccos(np.dot(mvec_ps_side, wvec_side)/(vnorm(mvec_ps_side) * vnorm(wvec_side))))
    assert (ideg_ps - jdeg_ps) > -1e-7

    # compute distances
    l1_ps = np.sqrt( np.linalg.norm(rvec)**2 + np.linalg.norm(mvec_ps_side)**2 - 2 * np.dot(rvec, mvec_ps_side) )
    l2_ps = np.sqrt( np.linalg.norm(mvec_ps_side)**2 + np.linalg.norm(wvec_side)**2 - 2 * np.dot(mvec_ps_side, wvec_side) )
    tps = l1_ps/cp + l2_ps/cs

    # Function to compute SP
    def func_sp(x):
        p, q = x
        if p<0 or q<0:
            return [np.inf, np.inf] 
        
        mvec = p*rvec + q*wvec_side
        svec = mvec - rvec
        cos_sp_j = np.dot(svec, wvec_side)/(vnorm(svec) * vnorm(wvec_side))
        cos_sp_i = np.dot(mvec, wvec_side)/(vnorm(mvec) * vnorm(wvec_side))
        sin_sp_j = np.sqrt(1-cos_sp_j**2)
        sin_sp_i = np.sqrt(1-cos_sp_i**2)
        f = np.zeros(2)
        f[0] = np.dot(mvec, wvec_side) / vnorm(wvec_side)**2 - 0.5
        f[1] = sin_sp_j/cs - sin_sp_i/cp
        return f

    B = fsolve(func_sp, [0.1, 0.1], full_output=True, xtol=1.49012e-08)
    print(B)
    # compute i and j
    p0_sp, q0_sp = B[0]
    mvec_sp_side = p0_sp*rvec + q0_ps*wvec_side
    svec_sp_side = mvec_sp_side - rvec

    jdeg_sp = np.rad2deg(np.arccos(np.dot(svec_sp_side, wvec_side)/(vnorm(svec_sp_side) * vnorm(wvec_side))))
    ideg_sp = np.rad2deg(np.arccos(np.dot(mvec_sp_side, wvec_side)/(vnorm(mvec_sp_side) * vnorm(wvec_side))))

    print(ideg_sp, jdeg_sp)
    assert (ideg_sp - jdeg_sp) > -1e-7

    # compute distances
    l1_sp = np.sqrt( np.linalg.norm(rvec)**2 + np.linalg.norm(mvec_sp_side)**2 - 2 * np.dot(rvec, mvec_sp_side) )
    l2_sp = np.sqrt( np.linalg.norm(mvec_sp_side)**2 + np.linalg.norm(wvec_side)**2 - 2 * np.dot(mvec_sp_side, wvec_side) )
    tsp = l1_sp/cs + l2_sp/cp

    return tps, tsp, (ideg_ps, jdeg_ps, ideg_sp, jdeg_sp)

def compute_VR(M0hat_try, st_event_trimmed, VR_sensors, weight_all, param):
    """
    compute VR associated with P and S window
    return VR with respect to only P, only S and P and S window
    """
    VR_comp_p = {} # dictionary to store VR associated with p window of each sensor
    VR_comp_s = {} # dictionary to store VR associated with s window of each sensor

    for sensorid in VR_sensors:
        sensorid_num = int(sensorid.split("OL")[-1])

        st_sensor = st_event_trimmed.select(station=sensorid).copy()

        tr_obs_Retrim_P =  st_sensor.select(location="Retrim_P", channel="OZ")[0]
        tr_syn_Retrim_P =  st_sensor.select(location="Retrim_P", channel="VZ")[0]
        tr_obs_Retrim_S =  st_sensor.select(location="Retrim_S", channel="OZ")[0]
        tr_syn_Retrim_S =  st_sensor.select(location="Retrim_S", channel="VZ")[0]

        # noramlize the signal with the maximum amplitude of obs
        
        
        # apply demean and taper
        
        tr_obs_Retrim_P.detrend(type="demean").taper(0.2)
        tr_syn_Retrim_P.detrend(type="demean").taper(0.2)
        tr_obs_Retrim_S.detrend(type="demean").taper(0.2)
        tr_syn_Retrim_S.detrend(type="demean").taper(0.2)
        
        # compute VR associated with P
        VRtmp1_P = np.linalg.norm(tr_obs_Retrim_P.data - M0hat_try * tr_syn_Retrim_P.data, axis=0) ** 2 # here the synthetic is scaled by M0hat_try
        VRtmp2_P = np.linalg.norm(tr_obs_Retrim_P.data, axis=0)**2
        # VRtmp2_P = np.linalg.norm(M0hat_try * tr_syn_Retrim_P.data, axis=0)**2
        VR_comp_p[sensorid] = weight_all[sensorid+"_P"] * (VRtmp1_P/VRtmp2_P)
        # VR_comp_p[sensorid] = weight_all[sensorid] * VRtmp1_P
        # Weighting with a factor of the absolute maximum amplitude of observation
        # weight_amp_p = param["weight_amp_factor"] * tr_obs_Retrim_P.max()
        # print(weight_amp_p)
        # VR_comp_p[sensorid] = weight_amp_p * (VRtmp1_P/VRtmp2_P)

        # compute VR associated with S
        VRtmp1_S = np.linalg.norm(tr_obs_Retrim_S.data - M0hat_try * tr_syn_Retrim_S.data, axis=0) ** 2 # here the synthetic is scaled by M0hat_try
        VRtmp2_S = np.linalg.norm(tr_obs_Retrim_S.data, axis=0)**2
        # VRtmp2_S = np.linalg.norm(M0hat_try * tr_syn_Retrim_S.data, axis=0)**2
        VR_comp_s[sensorid] = weight_all[sensorid+"_S"] * (VRtmp1_S/VRtmp2_S)
        # VR_comp_s[sensorid] = weight_all[sensorid] * VRtmp1_S
        # weight_amp_s = param["weight_amp_factor"] * tr_obs_Retrim_S.max()
        # VR_comp_s[sensorid] = weight_amp_s * (VRtmp1_S/VRtmp2_S)

    # print(VR_comp_p)
    # print(VR_comp_s)

    # sum up to compute VR
    VR_P  = 1.0 - sum(VR_comp_p.values())
    VR_S  = 1.0 - sum(VR_comp_s.values())
    VR_PS = 1.0 - param["alpha_p"] * sum(VR_comp_p.values()) - param["alpha_s"] * sum(VR_comp_s.values())
    
    return VR_P, VR_S, VR_PS, VR_comp_p, VR_comp_s



