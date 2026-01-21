#!/usr/bin/env python

# Functions to plot best fit traces
# 2022.05.09 Kurama Okubo

import os
import obspy
from obspy import read, Stream, Trace
from scipy import signal
import matplotlib.pyplot as plt
import glob
from glob import glob
import numpy as np
import pandas as pd
import datetime
from datetime import timedelta
from tqdm import tqdm
import pickle
import warnings
from obspy.core.utcdatetime import UTCDateTime    
import time
from scipy.spatial.distance import euclidean, correlation
from obspy.signal.cross_correlation import correlate, xcorr_max
from scipy.optimize import fsolve
from AEevents_Gridsearch_PS_func import *

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 12
os.environ['TZ'] = 'GMT' # change time zone to avoid confusion in unix_tvec conversion
UTCDateTime.DEFAULT_PRECISION = 8 # increase the time precision

def AE_GridSearch_plot_raw(st_event_trimmed, param, figsize=(11, 9)):
    """
    plot with shifting the trace using the P wave correlation
    """
    Nsensor = len(param["df_dist_sorted"])
    fig, axs = plt.subplots(Nsensor, 1, figsize=figsize, sharex=True)

    xlimit = [0, 0.12]
    # number of data point

    for i, sensorid in enumerate(param["df_dist_sorted"].keys()):

    # i = 0
    # sensorid = "OL15"
        sensorid_num = int(sensorid.split("OL")[-1])
        # select traces
        st_sensor = st_event_trimmed.select(station=sensorid).copy()
        # trim st_sensor within plotting range
        st_tmp = st_sensor[0].stats.starttime
        starttime_sensor = st_tmp
        endtime_sensor = st_tmp+xlimit[-1]*1e-3
        tr_obs_filt = st_sensor.select(location="Filt", channel="OZ")[0]
        tr_syn_filt = st_sensor.select(location="Filt", channel="VZ")[0]
        tr_syn_p_filt_shifted = st_sensor.select(location="FiltShifted_P", channel="VZ")[0]
        tr_obs_Retrim_P =  st_sensor.select(location="Retrim_P", channel="OZ")[0]
        tr_syn_Retrim_P =  st_sensor.select(location="Retrim_P", channel="VZ")[0]
        tr_obs_Retrim_S =  st_sensor.select(location="Retrim_S", channel="OZ")[0]
        tr_syn_Retrim_S =  st_sensor.select(location="Retrim_S", channel="VZ")[0]

        pretrigger = tr_obs_filt.stats.pretrigger #[ms]

        # 1. plot entire trace filter and filter_shifted
        axs[i].plot(tr_obs_filt.times()*1e3-pretrigger, tr_obs_filt.data*1e3, "k-", lw=1.2)
#         axs[i].plot(tr_syn_filt.times()*1e3-pretrigger, tr_syn_filt.data*1e3, "r")
        dt_shift_p = tr_syn_p_filt_shifted.stats.dt_shift
        print(dt_shift_p)
        # axs[i].plot((tr_syn_p_filt_shifted.times()+dt_shift_p + param["TR"]/2)*1e3-pretrigger, tr_syn_p_filt_shifted.data*1e3, "r-")
        axs[i].plot((tr_syn_p_filt_shifted.times()+dt_shift_p)*1e3-pretrigger, tr_syn_p_filt_shifted.data*1e3, "r-", lw=1.2)
        axs[i].set_xlim(xlimit)

       
        # set ylim
        fs = tr_obs_Retrim_S.stats.delta 
    #     npts_all = int(np.round( (xlimit[-1] + pretrigger)*1e-3/fs) )

        if "ylim" in param:
            axs[i].set_ylim([-param["ylim"][i], param["ylim"][i]])
            ylim_amp_filt = param["ylim"][i]
        else:
       #     ylim_amp_filt = 1.2*np.maximum(np.max(np.abs(tr_obs_filt.data[:npts_all]*1e3)), np.max(np.abs(tr_syn_p_filt_shifted.data[:npts_all]*1e3)))
            ylim_amp_filt = 3.0*np.maximum(np.max(np.abs(tr_obs_Retrim_S.data*1e3)), np.max(np.abs(tr_syn_Retrim_S.data*1e3)))
            axs[i].set_ylim([-ylim_amp_filt, ylim_amp_filt])

        event_id = param["datacase"].split("__")[1]
        axs[i].set_ylabel("{}\nvelocity [mm/s]".format(f"AS{sensorid[2:4]}"))
        # axs[0].set_title("event {} Bandpass filtered: {:.0f}-{:.0f}kHz".format(event_id, param["freqmin"]/1e3, param["freqmax"]/1e3))
        
        if "gougeevent_id" in param:
            axs[0].set_title(r"{} Mw={:4.1f} rake={:.1f}° TR={:.1f}μs Bandpass filtered: {:.0f}-{:.0f}kHz".format(
            param["gougeevent_id"], param["Mw"], param["rake"],  param["TR"]*1e6, param["freqmin"]/1e3, param["freqmax"]/1e3))
#             axs[0].set_title(" ")


        else:
            axs[0].set_title(r"Mw={:4.1f} rake={:.1f}° TR={:.1f}μs Bandpass filtered: {:.0f}-{:.0f}kHz".format(
            param["Mw"], param["rake"],  param["TR"]*1e6, param["freqmin"]/1e3, param["freqmax"]/1e3))

        # plot p and s window bar
    #     axs[i].plot(np.array([pwin_st, pwin_et])*1e3, 0.5*np.array([ylim_amp_filt, ylim_amp_filt]), "k-", lw=2)
    #     axs[i].plot(np.array([swin_st, swin_et])*1e3, 0.8*np.array([ylim_amp_filt, ylim_amp_filt]), "k-", lw=2)
    #     axs[i].text(pwin_st*1e3,  0.5*ylim_amp_filt, "P ", weight="normal", va="center", ha="right")
    #     axs[i].text(swin_st*1e3,  0.8*ylim_amp_filt, "S ", weight="normal", va="center", ha="right")

        # plot p and s arrival time
        p_arrival = tr_obs_filt.stats.p_arrival
        s_arrival = tr_obs_filt.stats.s_arrival
        axs[i].axvline(p_arrival+dt_shift_p*1e3, ls=":", c="k", lw=1)
        axs[i].axvline(s_arrival+dt_shift_p*1e3, ls=":", c="k", lw=1)
        # axs[i].axvline(p_arrival, ls=":", c="k", lw=1)
        # axs[i].axvline(s_arrival, ls=":", c="k", lw=1)
        
        if param["plot_reflection"]:
            # compute reflection
            p_direct, pp_side, pp_bottom, pp_top, ppp_side = compute_p_reflection(tr_obs_filt, param)
            tps, tsp, (ideg_ps, jdeg_ps, ideg_sp, jdeg_sp) = compute_ps_and_sp_side_reflection(tr_obs_filt, param)
            axs[i].axvline((pp_side+dt_shift_p)*1e3, c="g", ls="-", label="pp side", alpha=0.5)
            axs[i].axvline((pp_bottom+dt_shift_p)*1e3, c="m", ls="-", label="pp bottom", alpha=0.5)
            axs[i].axvline((pp_top+dt_shift_p)*1e3, c="c", ls="-", label="pp top", alpha=0.5)
            axs[i].axvline((ppp_side+dt_shift_p)*1e3, c="purple", ls="-", label="ppp side", alpha=0.5)
            axs[i].axvline((tps+dt_shift_p)*1e3, c="gray", ls="-", label="ps side", alpha=0.5)
            axs[i].axvline((tsp+dt_shift_p)*1e3, c="orange", ls="-", label="sp side", alpha=0.5)
        

        # plot p and s window bar
        # compute p and s window position (time after pretrigger)
        # pwin_st = tr_obs_Retrim_P.stats.pwin_start_tmp.microseconds * 1e-6  #[s]
        # pwin_et = tr_obs_Retrim_P.stats.pwin_end_tmp.microseconds * 1e-6  #[s]
        # swin_st = tr_obs_Retrim_S.stats.swin_start_tmp.microseconds * 1e-6  #[s]
        # swin_et = tr_obs_Retrim_S.stats.swin_end_tmp.microseconds * 1e-6  #[s]
        pwin_st = tr_obs_Retrim_P.stats.pwin_start_tmp  #[s]
        pwin_et = tr_obs_Retrim_P.stats.pwin_end_tmp  #[s]
        swin_st = tr_obs_Retrim_S.stats.swin_start_tmp  #[s]
        swin_et = tr_obs_Retrim_S.stats.swin_end_tmp  #[s]

        axs[i].plot(np.array([pwin_st+dt_shift_p, pwin_et+dt_shift_p])*1e3, 0.5*np.array([ylim_amp_filt, ylim_amp_filt]), "k-", lw=2)
        axs[i].plot(np.array([swin_st+dt_shift_p, swin_et+dt_shift_p])*1e3, 0.8*np.array([ylim_amp_filt, ylim_amp_filt]), "k-", lw=2)
        axs[i].text(pwin_st*1e3+dt_shift_p*1e3,  0.5*ylim_amp_filt, "P ", weight="normal", va="center", ha="right")
        axs[i].text(swin_st*1e3+dt_shift_p*1e3,  0.7*ylim_amp_filt, "S ", weight="normal", va="center", ha="right")

        # axs[i].get_yaxis().set_label_coords(-0.065,0.5)
        axs[i].get_yaxis().set_label_coords(-0.095,0.5)

    axs[-1].set_xlabel("Time [ms]")

    if param["plot_reflection"]:
        axs[0].legend(bbox_to_anchor=(1.04, 1.0), loc="upper left")
    else:
        axs[0].legend(["observation", "synthetic"], bbox_to_anchor=(1.04, 1.0), loc="upper left")
#         axs[0].legend(["observation"], bbox_to_anchor=(1.04, 1.0), loc="upper left")

def AE_GridSearch_plot_raw_noshift(st_event_trimmed, param):
    """
    plot only the non-shifted raw data
    """
    Nsensor = len(param["df_dist_sorted"])
    fig, axs = plt.subplots(Nsensor, 1, figsize=(9, 9), sharex=True)

    xlimit = [0, 0.1]
    # number of data point

    for i, sensorid in enumerate(param["df_dist_sorted"].keys()):

    # i = 0
    # sensorid = "OL15"
        sensorid_num = int(sensorid.split("OL")[-1])
        # select traces
        st_sensor = st_event_trimmed.select(station=sensorid).copy()
        # trim st_sensor within plotting range
        st_tmp = st_sensor[0].stats.starttime
        starttime_sensor = st_tmp
        endtime_sensor = st_tmp+xlimit[-1]*1e-3
        tr_obs_filt = st_sensor.select(location="Filt", channel="OZ")[0]
        tr_syn_filt = st_sensor.select(location="Filt", channel="VZ")[0]
        tr_syn_p_filt_shifted = st_sensor.select(location="FiltShifted_P", channel="VZ")[0]
        tr_obs_Retrim_P =  st_sensor.select(location="Retrim_P", channel="OZ")[0]
        tr_syn_Retrim_P =  st_sensor.select(location="Retrim_P", channel="VZ")[0]
        tr_obs_Retrim_S =  st_sensor.select(location="Retrim_S", channel="OZ")[0]
        tr_syn_Retrim_S =  st_sensor.select(location="Retrim_S", channel="VZ")[0]

        pretrigger = tr_obs_filt.stats.pretrigger #[ms]

        # 1. plot entire trace filter and filter_shifted
        axs[i].plot(tr_obs_filt.times()*1e3-pretrigger, tr_obs_filt.data*1e3, "k-", lw=1.2)
        axs[i].plot(tr_syn_filt.times()*1e3-pretrigger + param["TR"]*1e3/2, tr_syn_filt.data*1e3, "r", lw=1.2)
        dt_shift_p = tr_syn_p_filt_shifted.stats.dt_shift
        print(dt_shift_p)
        # axs[i].plot((tr_syn_p_filt_shifted.times()+dt_shift_p)*1e3-pretrigger, tr_syn_p_filt_shifted.data*1e3, "r-")
        axs[i].set_xlim(xlimit)

       
        # set ylim
        fs = tr_obs_Retrim_S.stats.delta 
    #     npts_all = int(np.round( (xlimit[-1] + pretrigger)*1e-3/fs) )
        if param["ylim"]:
            axs[i].set_ylim([-param["ylim"][i], param["ylim"][i]])
        else:
        #     ylim_amp_filt = 1.2*np.maximum(np.max(np.abs(tr_obs_filt.data[:npts_all]*1e3)), np.max(np.abs(tr_syn_p_filt_shifted.data[:npts_all]*1e3)))
            ylim_amp_filt = 2.0*np.maximum(np.max(np.abs(tr_obs_Retrim_S.data*1e3)), np.max(np.abs(tr_syn_Retrim_S.data*1e3)))
            axs[i].set_ylim([-ylim_amp_filt, ylim_amp_filt])

        event_id = param["datacase"].split("__")[1]
        axs[i].set_ylabel("{}\nvelocity [mm/s]".format(sensorid))
        axs[0].set_title(r"Mw={:4.1f} rake={:.1f}° TR={:.1f}μs Bandpass filtered: {:.0f}-{:.0f}kHz".format(
            param["Mw"], param["rake"],  param["TR"]*1e6, param["freqmin"]/1e3, param["freqmax"]/1e3))


        # plot p and s window bar
    #     axs[i].plot(np.array([pwin_st, pwin_et])*1e3, 0.5*np.array([ylim_amp_filt, ylim_amp_filt]), "k-", lw=2)
    #     axs[i].plot(np.array([swin_st, swin_et])*1e3, 0.8*np.array([ylim_amp_filt, ylim_amp_filt]), "k-", lw=2)
    #     axs[i].text(pwin_st*1e3,  0.5*ylim_amp_filt, "P ", weight="normal", va="center", ha="right")
    #     axs[i].text(swin_st*1e3,  0.8*ylim_amp_filt, "S ", weight="normal", va="center", ha="right")

        # plot p and s arrival time
        p_arrival = tr_obs_filt.stats.p_arrival
        s_arrival = tr_obs_filt.stats.s_arrival
        # axs[i].axvline(p_arrival+dt_shift_p*1e3, ls=":", c="k", lw=1)
        # axs[i].axvline(s_arrival+dt_shift_p*1e3, ls=":", c="k", lw=1)
        axs[i].axvline(p_arrival, ls=":", c="k", lw=1)
        axs[i].axvline(s_arrival, ls=":", c="k", lw=1)
        
        if param["plot_reflection"]:
            # compute reflection
            p_direct, pp_side, pp_bottom, pp_top, ppp_side = compute_p_reflection(tr_obs_filt, param)
            tps, tsp, (ideg_ps, jdeg_ps, ideg_sp, jdeg_sp) = compute_ps_and_sp_side_reflection(tr_obs_filt, param)
            axs[i].axvline((pp_side)*1e3, c="g", ls="-", label="pp side", alpha=0.5)
            axs[i].axvline((pp_bottom)*1e3, c="m", ls="-", label="pp bottom", alpha=0.5)
            axs[i].axvline((pp_top)*1e3, c="c", ls="-", label="pp top", alpha=0.5)
            axs[i].axvline((ppp_side)*1e3, c="purple", ls="-", label="ppp side", alpha=0.5)
            axs[i].axvline((tps)*1e3, c="gray", ls="-", label="ps side", alpha=0.5)
            axs[i].axvline((tsp)*1e3, c="orange", ls="-", label="sp side", alpha=0.5)
        

        # plot p and s window bar
        # compute p and s window position (time after pretrigger)
        # pwin_st = tr_obs_Retrim_P.stats.pwin_start_tmp.microseconds * 1e-6  #[s]
        # pwin_et = tr_obs_Retrim_P.stats.pwin_end_tmp.microseconds * 1e-6  #[s]
        # swin_st = tr_obs_Retrim_S.stats.swin_start_tmp.microseconds * 1e-6  #[s]
        # swin_et = tr_obs_Retrim_S.stats.swin_end_tmp.microseconds * 1e-6  #[s]
        # pwin_st = tr_obs_Retrim_P.stats.pwin_start_tmp  #[s]
        # pwin_et = tr_obs_Retrim_P.stats.pwin_end_tmp  #[s]
        # swin_st = tr_obs_Retrim_S.stats.swin_start_tmp  #[s]
        # swin_et = tr_obs_Retrim_S.stats.swin_end_tmp  #[s]

        # axs[i].plot(np.array([pwin_st, pwin_et])*1e3, 0.5*np.array([ylim_amp_filt, ylim_amp_filt]), "k-", lw=2)
        # axs[i].plot(np.array([swin_st, swin_et])*1e3, 0.8*np.array([ylim_amp_filt, ylim_amp_filt]), "k-", lw=2)
        # axs[i].text(pwin_st*1e3,  0.5*ylim_amp_filt, "P ", weight="normal", va="center", ha="right")
        # axs[i].text(swin_st*1e3,  0.8*ylim_amp_filt, "S ", weight="normal", va="center", ha="right")

        # axs[i].get_yaxis().set_label_coords(-0.08,0.5)
        axs[i].get_yaxis().set_label_coords(-0.095,0.5)

    # fig.text(0.04, 0.5, 'Velocity [mm/s]', va='center', rotation='vertical')

    axs[-1].set_xlabel("Time [ms]")

    if param["plot_reflection"]:
        axs[0].legend(bbox_to_anchor=(1.01, 1.0), loc="upper left")
    else:
        axs[0].legend(["observation", "synthetic"], bbox_to_anchor=(1.01, 1.0), loc="upper left")

def AE_GridSearch_plot_reflection(st_event_trimmed, param, i = 2):
    """
    plot only the non-shifted raw data
    """
    Nsensor = len(param["df_dist_sorted"])
    # fig, axs = plt.subplots(Nsensor, 1, figsize=(9, 9), sharex=True)
    fig, ax = plt.subplots(1, 1, figsize=(8, 3))

    xlimit = [0, 0.1]
    # number of data point
    sensorids = list(param["df_dist_sorted"].keys())
    # for i, sensorid in enumerate(param["df_dist_sorted"].keys()):

    sensorid = sensorids[i]
    # sensorid = "OL15"
    sensorid_num = int(sensorid.split("OL")[-1])
    # select traces
    st_sensor = st_event_trimmed.select(station=sensorid).copy()
    # trim st_sensor within plotting range
    st_tmp = st_sensor[0].stats.starttime
    starttime_sensor = st_tmp
    endtime_sensor = st_tmp+xlimit[-1]*1e-3
    tr_obs_filt = st_sensor.select(location="Filt", channel="OZ")[0]
    tr_syn_filt = st_sensor.select(location="Filt", channel="VZ")[0]
    tr_syn_p_filt_shifted = st_sensor.select(location="FiltShifted_P", channel="VZ")[0]
    tr_obs_Retrim_P =  st_sensor.select(location="Retrim_P", channel="OZ")[0]
    tr_syn_Retrim_P =  st_sensor.select(location="Retrim_P", channel="VZ")[0]
    tr_obs_Retrim_S =  st_sensor.select(location="Retrim_S", channel="OZ")[0]
    tr_syn_Retrim_S =  st_sensor.select(location="Retrim_S", channel="VZ")[0]

    pretrigger = tr_obs_filt.stats.pretrigger #[ms]

    # 1. plot entire trace filter and filter_shifted
    # axs[i].plot(tr_obs_filt.times()*1e3-pretrigger, tr_obs_filt.data*1e3, "k-")
    ax.plot(tr_syn_filt.times()*1e3-pretrigger + param["TR"]*1e3/2, tr_syn_filt.data*1e3, "r", label="synthetic")
    dt_shift_p = tr_syn_p_filt_shifted.stats.dt_shift
    print(dt_shift_p)
    # ax.plot((tr_syn_p_filt_shifted.times()+dt_shift_p)*1e3-pretrigger, tr_syn_p_filt_shifted.data*1e3, "r-")
    ax.set_xlim(xlimit)

   
    # set ylim
    fs = tr_obs_Retrim_S.stats.delta 
#     npts_all = int(np.round( (xlimit[-1] + pretrigger)*1e-3/fs) )
#     ylim_amp_filt = 1.2*np.maximum(np.max(np.abs(tr_obs_filt.data[:npts_all]*1e3)), np.max(np.abs(tr_syn_p_filt_shifted.data[:npts_all]*1e3)))
    ylim_amp_filt = 3.0*np.maximum(np.max(np.abs(tr_obs_Retrim_S.data*1e3)), np.max(np.abs(tr_syn_Retrim_S.data*1e3)))
    ax.set_ylim([-ylim_amp_filt, ylim_amp_filt])

    event_id = param["datacase"].split("__")[1]
    # ax.set_ylabel("{}\nvelocity [mm/s]".format(sensorid))
    ax.set_title("{} Bandpass filtered: {:.0f}-{:.0f}kHz".format(sensorid, param["freqmin"]/1e3, param["freqmax"]/1e3))


    # plot p and s window bar
#     ax.plot(np.array([pwin_st, pwin_et])*1e3, 0.5*np.array([ylim_amp_filt, ylim_amp_filt]), "k-", lw=2)
#     ax.plot(np.array([swin_st, swin_et])*1e3, 0.8*np.array([ylim_amp_filt, ylim_amp_filt]), "k-", lw=2)
#     ax.text(pwin_st*1e3,  0.5*ylim_amp_filt, "P ", weight="normal", va="center", ha="right")
#     ax.text(swin_st*1e3,  0.8*ylim_amp_filt, "S ", weight="normal", va="center", ha="right")

    # plot p and s arrival time
    p_arrival = tr_obs_filt.stats.p_arrival
    s_arrival = tr_obs_filt.stats.s_arrival
    # ax.axvline(p_arrival+dt_shift_p*1e3, ls=":", c="k", lw=1)
    # ax.axvline(s_arrival+dt_shift_p*1e3, ls=":", c="k", lw=1)
    ax.axvline(p_arrival, ls="-", c="k", lw=1.2, zorder = -1, label="p, s")
    ax.axvline(s_arrival, ls="-", c="k", lw=1.2, zorder = -1)
    
    if param["plot_reflection"]:
        # compute reflection
        p_direct, pp_side, pp_bottom, pp_top, ppp_side = compute_p_reflection(tr_obs_filt, param)
        tps, tsp, (ideg_ps, jdeg_ps, ideg_sp, jdeg_sp) = compute_ps_and_sp_side_reflection(tr_obs_filt, param)
        ax.axvline((pp_side)*1e3, c="g", ls="-", label="pp side", alpha=0.5)
        ax.axvline((pp_bottom)*1e3, c="m", ls="-", label="pp bottom", alpha=0.5)
        ax.axvline((pp_top)*1e3, c="c", ls="-", label="pp top", alpha=0.5)
        ax.axvline((ppp_side)*1e3, c="purple", ls="-", label="ppp side", alpha=0.5)
        ax.axvline((tps)*1e3, c="gray", ls="-", label="ps side", alpha=0.5)
        ax.axvline((tsp)*1e3, c="orange", ls="-", label="sp side", alpha=0.5)
    

    # plot p and s window bar
    # compute p and s window position (time after pretrigger)
    # pwin_st = tr_obs_Retrim_P.stats.pwin_start_tmp.microseconds * 1e-6  #[s]
    # pwin_et = tr_obs_Retrim_P.stats.pwin_end_tmp.microseconds * 1e-6  #[s]
    # swin_st = tr_obs_Retrim_S.stats.swin_start_tmp.microseconds * 1e-6  #[s]
    # swin_et = tr_obs_Retrim_S.stats.swin_end_tmp.microseconds * 1e-6  #[s]
    # pwin_st = tr_obs_Retrim_P.stats.pwin_start_tmp  #[s]
    # pwin_et = tr_obs_Retrim_P.stats.pwin_end_tmp  #[s]
    # swin_st = tr_obs_Retrim_S.stats.swin_start_tmp  #[s]
    # swin_et = tr_obs_Retrim_S.stats.swin_end_tmp  #[s]

    # axs[i].plot(np.array([pwin_st, pwin_et])*1e3, 0.5*np.array([ylim_amp_filt, ylim_amp_filt]), "k-", lw=2)
    # axs[i].plot(np.array([swin_st, swin_et])*1e3, 0.8*np.array([ylim_amp_filt, ylim_amp_filt]), "k-", lw=2)
    # axs[i].text(pwin_st*1e3,  0.5*ylim_amp_filt, "P ", weight="normal", va="center", ha="right")
    # axs[i].text(swin_st*1e3,  0.8*ylim_amp_filt, "S ", weight="normal", va="center", ha="right")

    # axs[i].get_yaxis().set_label_coords(-0.08,0.5)

    fig.text(0.04, 0.5, 'Velocity [mm/s]', va='center', rotation='vertical')

    ax.set_xlabel("Time [ms]")

    if param["plot_reflection"]:
        ax.legend(bbox_to_anchor=(1.04, 1.0), loc="upper left")
    else:
        ax.legend(["observation", "synthetic"], bbox_to_anchor=(1.04, 1.0), loc="upper left")

def AE_GridSearch_plot_entirewaveform_withP(st_event_trimmed, param):
    # plot comparison of waveform
    # NOTE: tr_filt includes pretigger, while tr_retrimmed is retrimmed without pretigger. Thus, correct pretrigger only on the filter (entire) trace.
    # pwin_st, pwin_et does not contain the pretrigger. 

    Nsensor = len(param["df_dist_sorted"])

    fig, axs = plt.subplots(Nsensor, 2, figsize=(10.5, 10.8), gridspec_kw={'width_ratios': [5, 1]})
#     fig, axs = plt.subplots(Nsensor, 2, figsize=(12, 10.8), gridspec_kw={'width_ratios': [5, 1]})

    xlimit = [0, 0.12]
    # number of data point

    for i, sensorid in enumerate(param["df_dist_sorted"].keys()):

        # i = 0
        # sensorid = "OL15"
        sensorid_num = int(sensorid.split("OL")[-1])
        # select traces
        st_sensor = st_event_trimmed.select(station=sensorid).copy()
        # trim st_sensor within plotting range
        st_tmp = st_sensor[0].stats.starttime
        starttime_sensor = st_tmp
        endtime_sensor = st_tmp+xlimit[-1]*1e-3
        tr_obs_filt = st_sensor.select(location="Filt", channel="OZ")[0]
        tr_syn_filt = st_sensor.select(location="Filt", channel="VZ")[0]
        tr_syn_p_filt_shifted = st_sensor.select(location="FiltShifted_P", channel="VZ")[0]
        tr_obs_Retrim_P =  st_sensor.select(location="Retrim_P", channel="OZ")[0]
        tr_syn_Retrim_P =  st_sensor.select(location="Retrim_P", channel="VZ")[0]
        tr_obs_Retrim_S =  st_sensor.select(location="Retrim_S", channel="OZ")[0]
        tr_syn_Retrim_S =  st_sensor.select(location="Retrim_S", channel="VZ")[0]

        pretrigger = tr_obs_filt.stats.pretrigger #[ms]

        # 1. plot entire trace filter and filter_shifted
        # axs[i, 0].plot(tr_obs_filt.times()*1e3-pretrigger + param["TR"]*1e3/2, tr_obs_filt.data*1e3, "k-")

        # NOTE: the aperture correction on the entire waveform screwed up the amplitude on the reflected P waves, which would have different incident angle.
        # if param['aperturecorrection']:
        #     # divided the amplitude by the factor of aperture effect
        #     df_incidentangle_sensor = param["df_incidentangle"].loc[f"{sensorid}__{param['datacase']}"]
        #     incidentangle = df_incidentangle_sensor["incidentangle"]
        #     beta_coef_p = float(incidentangle_scalingfactor_analytic(param["cp"], np.deg2rad(incidentangle), param["TR"], param["R_sensor"]))
        #     # print(f"aperture correct: {sensorid} incident angle={df_incidentangle_sensor['incidentangle']:4f}deg. beta_p={beta_coef_p:4f}")   
        # else:
        #     beta_coef_p = 1.0

        axs[i, 0].plot(tr_obs_filt.times()*1e3-pretrigger, (tr_obs_filt.data)*1e3, "k-", lw=1.2)
        # axs[i, 0].plot(tr_syn_filt.times()*1e3-pretrigger, tr_syn_filt.data*1e3, "g:")
        dt_shift_p = tr_syn_Retrim_P.stats.dt_shift
        dt_shift_s = tr_syn_Retrim_S.stats.dt_shift
        # print(dt_shift)
        # axs[i].plot((tr_syn_p_filt_shifted.times()+dt_shift_p)*1e3-pretrigger, tr_syn_p_filt_shifted.data*1e3, "r-", lw=1.2)

        axs[i, 0].plot((tr_syn_p_filt_shifted.times()+dt_shift_p)*1e3-pretrigger, tr_syn_p_filt_shifted.data*1e3, "r-", lw=1.2)
        axs[i, 0].set_xlim(xlimit)        

        # compute p and s window position after shifted(time after pretrigger)
        # pwin_st = tr_obs_Retrim_P.stats.pwin_start.microseconds * 1e-6  #[ms]
        # pwin_et = tr_obs_Retrim_P.stats.pwin_end.microseconds * 1e-6  #[ms]
        # swin_st = tr_obs_Retrim_S.stats.swin_start.microseconds * 1e-6  #[ms]
        # swin_et = tr_obs_Retrim_S.stats.swin_end.microseconds * 1e-6  #[ms]

        pwin_st = tr_obs_Retrim_P.stats.pwin_start  #[s]
        pwin_et = tr_obs_Retrim_P.stats.pwin_end  #[s]
        swin_st = tr_obs_Retrim_S.stats.swin_start  #[s]
        swin_et = tr_obs_Retrim_S.stats.swin_end  #[s]

        # superimpose the shifted synthetic traces
        # axs[i, 0].plot((tr_syn_Retrim_P.times()+pwin_st+param["TR"]/2)*1e3, tr_syn_Retrim_P.data*1e3, "r-", lw=2)
        # axs[i, 0].plot((tr_syn_Retrim_S.times()+swin_st+param["TR"]/2)*1e3, tr_syn_Retrim_S.data*1e3, "r-", lw=2)
        # axs[i, 0].plot((tr_syn_Retrim_P.times()+pwin_st)*1e3, tr_syn_Retrim_P.data*1e3, "r-", lw=2)
        # axs[i, 0].plot((tr_syn_Retrim_S.times()+swin_st)*1e3, tr_syn_Retrim_S.data*1e3, "r-", lw=2)

        # 2. plot P window

        # threshold with correlation coefficient
        ls_P = "-"
        if "cc_threshold" in param:
            if tr_syn_Retrim_P.stats.cc.max() < param["cc_threshold"]:
                ls_P = "--"

        axs[i, 1].plot((tr_obs_Retrim_P.times()+pwin_st)*1e3, tr_obs_Retrim_P.detrend(type="demean").taper(0.05).data*1e3, "k-")
        axs[i, 1].plot((tr_syn_Retrim_P.times()+pwin_st)*1e3, tr_syn_Retrim_P.detrend(type="demean").taper(0.05).data*1e3, "r", ls=ls_P)

        # Plotting the P waveform without the aperture correction
        # re-multiply the beta_coef_p as the retrimmed observation is already scaled
        if "aperturecorrection" in param:
            # divided the amplitude by the factor of aperture effect
            df_incidentangle_sensor = param["df_incidentangle"].loc[f"{sensorid}__{param['datacase']}"]
            incidentangle = df_incidentangle_sensor["incidentangle"]
            beta_coef_p = float(incidentangle_scalingfactor_analytic(param["cp"], np.deg2rad(incidentangle), param["TR"], param["R_sensor"]))
            # print(f"aperture correct: {sensorid} incident angle={df_incidentangle_sensor['incidentangle']:4f}deg. beta_p={beta_coef_p:4f}")   
            axs[i, 1].plot((tr_obs_Retrim_P.times()+pwin_st)*1e3, tr_obs_Retrim_P.detrend(type="demean").taper(0.05).data*beta_coef_p*1e3, "b:")
        else:
            incidentangle = ""

        # # 3. plot S window
        # ls_S = "-"
        # if "cc_threshold" in param:
        #     if tr_syn_Retrim_S.stats.cc.max() < param["cc_threshold"]:
        #         ls_S = "--"

        # axs[i, 2].plot((tr_obs_Retrim_S.times()+swin_st)*1e3, tr_obs_Retrim_S.detrend(type="demean").taper(0.05).data*1e3, "k-")
        # axs[i, 2].plot((tr_syn_Retrim_S.times()+swin_st)*1e3, tr_syn_Retrim_S.detrend(type="demean").taper(0.05).data*1e3, "r", ls=ls_S)

        # set ylim
        # fs = tr_obs_Retrim_S.stats.delta 

    #     npts_all = int(np.round( (xlimit[-1] + pretrigger)*1e-3/fs) )
    #     ylim_amp_filt = 1.2*np.maximum(np.max(np.abs(tr_obs_filt.data[:npts_all]*1e3)), np.max(np.abs(tr_syn_p_filt_shifted.data[:npts_all]*1e3)))
        # ylim_amp_filt = 2.0*np.maximum(np.max(np.abs(tr_obs_Retrim_S.data*1e3)), np.max(np.abs(tr_syn_Retrim_S.data*1e3)))
        ylim_amp_filt_P = 4.0*np.maximum(np.max(np.abs(tr_obs_Retrim_P.data*1e3)), np.max(np.abs(tr_syn_Retrim_P.data*1e3)))
        
        if "ylim" in param:
            ylim_amp_filt_entirewaveform = param["ylim"][i]
            axs[i, 0].set_ylim([-ylim_amp_filt_entirewaveform, ylim_amp_filt_entirewaveform])
        else:
            axs[i, 0].set_ylim([-ylim_amp_filt_P, ylim_amp_filt_P])
            
        axs[i, 0].tick_params(axis='x', pad=7)

        axs[i, 1].set_ylim(np.array([-ylim_amp_filt_P, ylim_amp_filt_P])*0.5)

        # ylim_amp_filt_S = 1.5*np.maximum(np.max(np.abs(tr_obs_Retrim_S.data*1e3)), np.max(np.abs(tr_syn_Retrim_S.data*1e3)))
        # axs[i, 2].set_ylim(np.array([-ylim_amp_filt_S, ylim_amp_filt_S]))

        axs[i, 0].set_ylabel("AS{}\nvelocity [mm/s]".format(sensorid[2:]))
        event_id = param["datacase"].split("__")[1]
        if "M0" in param:
            title_x = [event_id, st_event_trimmed.select(channel="OZ")[0].stats.origintime, param["M0"],  param["TR"]*1e6, param["freqmin"]/1e3, param["freqmax"]/1e3, "M$_{\mathrm{0}}$", "T$_{\mathrm{R}}$"]
            titlestr = 'Event {0}: Origin time: {1:.4f}s\n{6}: {2:.2f}Nm {7}: {3:.1f}μs Bandpass filtered: {4:.0f}-{5:.0f}kHz'.format(*title_x)
        else:
            title_x = [event_id, st_event_trimmed.select(channel="OZ")[0].stats.origintime, param["freqmin"]/1e3, param["freqmax"]/1e3]
            titlestr = 'Event {0}: Origin time: {1:.4f}s Bandpass filtered: {2:.0f}-{3:.0f}kHz'.format(*title_x)  
            
        axs[0, 0].set_title(titlestr, fontsize=13)
        axs[0, 1].set_title("P wave")
        # axs[0, 2].set_title("S wave")

        # plot p and s window bar

        # axs[i, 0].plot(np.array([swin_st, swin_et])*1e3, 0.8*np.array([ylim_amp_filt_P, ylim_amp_filt_P]), "k-", lw=2)
        # axs[i, 0].text(pwin_st*1e3,  0.5*ylim_amp_filt_P, "P ", weight="normal", va="center", ha="right")
        # axs[i, 0].text(swin_st*1e3,  0.8*ylim_amp_filt_P, "S ", weight="normal", va="center", ha="right")

        # plot p and s arrival time
        p_arrival = tr_obs_filt.stats.p_arrival
        s_arrival = tr_obs_filt.stats.s_arrival
        axs[i, 0].axvline(p_arrival + dt_shift_p*1e3 , ls=":", c="k", lw=1)
        axs[i, 0].axvline(s_arrival + dt_shift_s*1e3 , ls=":", c="k", lw=1)
        if "ylim" in param:
            axs[i, 0].plot(np.array([pwin_st, pwin_et])*1e3, 0.5*np.array([ylim_amp_filt_entirewaveform, ylim_amp_filt_entirewaveform]), "k-", lw=2)
            axs[i, 0].text(p_arrival + dt_shift_p*1e3, 0.7*ylim_amp_filt_entirewaveform, " P",)
            axs[i, 0].text(s_arrival + dt_shift_s*1e3, 0.7*ylim_amp_filt_entirewaveform, " S",)

        else:
            axs[i, 0].plot(np.array([pwin_st, pwin_et])*1e3, 0.5*np.array([ylim_amp_filt_P, ylim_amp_filt_P]), "k-", lw=2)
            axs[i, 0].text(p_arrival + dt_shift_p*1e3, 0.7*ylim_amp_filt_P, " P",)
            axs[i, 0].text(s_arrival + dt_shift_s*1e3, 0.7*ylim_amp_filt_P, " S",)
 
        axs[i, 0].get_yaxis().set_label_coords(-0.0612,0.5)

        # Annotate the source distance and incident angle
        annot_x = [param['df_dist_sorted'][sensorid], incidentangle]
        if "ylim" in param:
            axs[i, 0].text(0.0895, -0.88*ylim_amp_filt_entirewaveform, "source distance: {0:.1f}mm\nincident angle:    {1:.1f}°".format(*annot_x), fontsize=10.5)
        else:
            axs[i, 0].text(0.0895, -0.88*ylim_amp_filt_P, "source distance: {0:.1f}mm\nincident angle:    {1:.1f}°".format(*annot_x), fontsize=10.5)

    axs[-1, 0].set_xlabel("Time [ms]")

    fig.tight_layout()

    
# def AE_GridSearch_plot_shifted(st_event_trimmed, param):
#     # plot comparison of waveform
#     # NOTE: tr_filt includes pretigger, while tr_retrimmed is retrimmed without pretigger. Thus, correct pretrigger only on the filter (entire) trace.
#     # pwin_st, pwin_et does not contain the pretrigger. 

#     Nsensor = len(param["df_dist_sorted"])

#     fig, axs = plt.subplots(Nsensor, 3, figsize=(12, 9), gridspec_kw={'width_ratios': [4, 1, 1]})

#     xlimit = [0, 0.12]
#     # number of data point

#     for i, sensorid in enumerate(param["df_dist_sorted"].keys()):

#         # i = 0
#         # sensorid = "OL15"
#         sensorid_num = int(sensorid.split("OL")[-1])
#         # select traces
#         st_sensor = st_event_trimmed.select(station=sensorid).copy()
#         # trim st_sensor within plotting range
#         st_tmp = st_sensor[0].stats.starttime
#         starttime_sensor = st_tmp
#         endtime_sensor = st_tmp+xlimit[-1]*1e-3
#         tr_obs_filt = st_sensor.select(location="Filt", channel="OZ")[0]
#         tr_syn_filt = st_sensor.select(location="Filt", channel="VZ")[0]
#         tr_syn_p_filt_shifted = st_sensor.select(location="FiltShifted_P", channel="VZ")[0]
#         tr_obs_Retrim_P =  st_sensor.select(location="Retrim_P", channel="OZ")[0]
#         tr_syn_Retrim_P =  st_sensor.select(location="Retrim_P", channel="VZ")[0]
#         tr_obs_Retrim_S =  st_sensor.select(location="Retrim_S", channel="OZ")[0]
#         tr_syn_Retrim_S =  st_sensor.select(location="Retrim_S", channel="VZ")[0]

#         pretrigger = tr_obs_filt.stats.pretrigger #[ms]

#         # 1. plot entire trace filter and filter_shifted
#         # axs[i, 0].plot(tr_obs_filt.times()*1e3-pretrigger + param["TR"]*1e3/2, tr_obs_filt.data*1e3, "k-")
#         axs[i, 0].plot(tr_obs_filt.times()*1e3-pretrigger, tr_obs_filt.data*1e3, "k-")
#         # axs[i, 0].plot(tr_syn_filt.times()*1e3-pretrigger, tr_syn_filt.data*1e3, "g:")
#         dt_shift_p = tr_syn_Retrim_P.stats.dt_shift
#         dt_shift_s = tr_syn_Retrim_S.stats.dt_shift
#         # print(dt_shift)
#         # axs[i, 0].plot((tr_syn_p_filt_shifted.times()+dt_shift_p + param["TR"]/2)*1e3-pretrigger, tr_syn_p_filt_shifted.data*1e3, "g:")
#         axs[i, 0].set_xlim(xlimit)


#         # compute p and s window position after shifted(time after pretrigger)
#         # pwin_st = tr_obs_Retrim_P.stats.pwin_start.microseconds * 1e-6  #[ms]
#         # pwin_et = tr_obs_Retrim_P.stats.pwin_end.microseconds * 1e-6  #[ms]
#         # swin_st = tr_obs_Retrim_S.stats.swin_start.microseconds * 1e-6  #[ms]
#         # swin_et = tr_obs_Retrim_S.stats.swin_end.microseconds * 1e-6  #[ms]

#         pwin_st = tr_obs_Retrim_P.stats.pwin_start  #[s]
#         pwin_et = tr_obs_Retrim_P.stats.pwin_end  #[s]
#         swin_st = tr_obs_Retrim_S.stats.swin_start  #[s]
#         swin_et = tr_obs_Retrim_S.stats.swin_end  #[s]

#         # superimpose the shifted synthetic traces
#         # axs[i, 0].plot((tr_syn_Retrim_P.times()+pwin_st+param["TR"]/2)*1e3, tr_syn_Retrim_P.data*1e3, "r-", lw=2)
#         # axs[i, 0].plot((tr_syn_Retrim_S.times()+swin_st+param["TR"]/2)*1e3, tr_syn_Retrim_S.data*1e3, "r-", lw=2)
#         axs[i, 0].plot((tr_syn_Retrim_P.times()+pwin_st)*1e3, tr_syn_Retrim_P.data*1e3, "r-", lw=2)
#         axs[i, 0].plot((tr_syn_Retrim_S.times()+swin_st)*1e3, tr_syn_Retrim_S.data*1e3, "r-", lw=2)

#         # 2. plot P window

#         # threshold with correlation coefficient
#         ls_P = "-"
#         if "cc_threshold" in param:
#             if tr_syn_Retrim_P.stats.cc.max() < param["cc_threshold"]:
#                 ls_P = "--"

#         axs[i, 1].plot((tr_obs_Retrim_P.times()+pwin_st)*1e3, tr_obs_Retrim_P.detrend(type="demean").taper(0.2).data*1e3, "k-")
#         axs[i, 1].plot((tr_syn_Retrim_P.times()+pwin_st)*1e3, tr_syn_Retrim_P.detrend(type="demean").taper(0.2).data*1e3, "r", ls=ls_P)

#         # 3. plot S window
#         ls_S = "-"
#         if "cc_threshold" in param:
#             if tr_syn_Retrim_S.stats.cc.max() < param["cc_threshold"]:
#                 ls_S = "--"

#         axs[i, 2].plot((tr_obs_Retrim_S.times()+swin_st)*1e3, tr_obs_Retrim_S.detrend(type="demean").taper(0.2).data*1e3, "k-")
#         axs[i, 2].plot((tr_syn_Retrim_S.times()+swin_st)*1e3, tr_syn_Retrim_S.detrend(type="demean").taper(0.2).data*1e3, "r", ls=ls_S)

#         # set ylim
#         fs = tr_obs_Retrim_S.stats.delta 

#     #     npts_all = int(np.round( (xlimit[-1] + pretrigger)*1e-3/fs) )
#     #     ylim_amp_filt = 1.2*np.maximum(np.max(np.abs(tr_obs_filt.data[:npts_all]*1e3)), np.max(np.abs(tr_syn_p_filt_shifted.data[:npts_all]*1e3)))
#         # ylim_amp_filt = 2.0*np.maximum(np.max(np.abs(tr_obs_Retrim_S.data*1e3)), np.max(np.abs(tr_syn_Retrim_S.data*1e3)))
#         ylim_amp_filt_P = 4.0*np.maximum(np.max(np.abs(tr_obs_Retrim_P.data*1e3)), np.max(np.abs(tr_syn_Retrim_P.data*1e3)))
#         axs[i, 0].set_ylim([-ylim_amp_filt_P, ylim_amp_filt_P])
#         axs[i, 1].set_ylim(np.array([-ylim_amp_filt_P, ylim_amp_filt_P])*0.5)

#         ylim_amp_filt_S = 1.5*np.maximum(np.max(np.abs(tr_obs_Retrim_S.data*1e3)), np.max(np.abs(tr_syn_Retrim_S.data*1e3)))
#         axs[i, 2].set_ylim(np.array([-ylim_amp_filt_S, ylim_amp_filt_S]))

#         axs[i, 0].set_ylabel("{}\nvelocity [mm/s]".format(sensorid))
#         event_id = param["datacase"].split("__")[1]
#         axs[0, 0].set_title("event {} Bandpass filtered: {:.0f}-{:.0f}kHz".format(event_id, param["freqmin"]/1e3, param["freqmax"]/1e3))
#         axs[0, 1].set_title("P wave")
#         axs[0, 2].set_title("S wave")

#         # plot p and s window bar

#         axs[i, 0].plot(np.array([pwin_st, pwin_et])*1e3, 0.5*np.array([ylim_amp_filt_P, ylim_amp_filt_P]), "k-", lw=2)
#         axs[i, 0].plot(np.array([swin_st, swin_et])*1e3, 0.8*np.array([ylim_amp_filt_P, ylim_amp_filt_P]), "k-", lw=2)
#         axs[i, 0].text(pwin_st*1e3,  0.5*ylim_amp_filt_P, "P ", weight="normal", va="center", ha="right")
#         axs[i, 0].text(swin_st*1e3,  0.8*ylim_amp_filt_P, "S ", weight="normal", va="center", ha="right")

#         # plot p and s arrival time
#         p_arrival = tr_obs_filt.stats.p_arrival
#         s_arrival = tr_obs_filt.stats.s_arrival
#         # axs[i, 0].axvline(p_arrival + dt_shift_p*1e3 , ls=":", c="k", lw=0.5)
#         # axs[i, 0].axvline(s_arrival + dt_shift_s*1e3 , ls=":", c="k", lw=0.5)

#         axs[i, 0].get_yaxis().set_label_coords(-0.095,0.5)


#     axs[-1, 0].set_xlabel("Time [ms]")

#     fig.tight_layout()


# define function to plot colorbar
# reference: https://pythonmatplotlibtips.blogspot.com/2019/07/draw-two-axis-to-one-colorbar.html
def plot_twincbar(cb, vmax):
    pos = cb.ax.get_position()
    ax1 = cb.ax
    ax1.set_aspect('auto')
    # create a second axis and specify ticks based on the relation between the first axis and second aces
    ax2 = ax1.twiny()
    newlabel = [0.8, 0.85, 0.9, 0.95, 1.0] # labels of the ticklabels: the position in the new axis
    normVR2VR= lambda x: x*vmax # convert function
    newpos   = [normVR2VR(x) for x in newlabel]  # position of the ticklabels in the old axis
    ax2.set_xlim([newpos[0], newpos[-1]])
    ax2.set_xticks(newpos)
    ax2.set_xticklabels(newlabel)

    # # resize the colorbar
    pos.y0 += 0.05
    pos.y1 -= 0.022
    pos.x0 += 0.15
    pos.x1 -= 0.15
    # # arrange and adjust the position of each axis, ticks, and ticklabels
    ax1.set_position(pos)
    ax2.set_position(pos)
    ax1.xaxis.set_ticks_position('top') # set the position of the first axis to right
    ax1.xaxis.set_label_position('top') # set the position of the fitst axis to right
    ax1.set_xlabel('VR')
    ax2.xaxis.set_ticks_position('bottom') # set the position of the second axis to right
    ax2.xaxis.set_label_position('bottom') # set the position of the second axis to right
    # ax2.spines['left'].set_position(('outward', 50)) # adjust the position of the second axis
    ax2.set_xlabel(r'$\overline{VR}$')
    
# define function to plot colorbar
def plot_twincbar_M0(cb, vmax):
    pos = cb.ax.get_position()
    ax1 = cb.ax
    ax1.set_aspect('auto')
    # create a second axis and specify ticks based on the relation between the first axis and second aces
    # ref: https://pythonmatplotlibtips.blogspot.com/2018/01/add-second-x-axis-below-first-x-axis-python-matplotlib-pyplot.html
    ax2 = ax1.twiny()
    newlabel = [0.8, 0.85, 0.9, 0.95, 1.0] # labels of the ticklabels: the position in the new axis
    normVR2VR= lambda x: x*vmax # convert function: from Kelvin to Degree Celsius
    newpos   = [normVR2VR(x) for x in newlabel]  # position of the ticklabels in the old axis
    ax2.set_xlim([newpos[0], newpos[-1]])
    ax2.set_xticks(newpos)
    ax2.set_xticklabels(newlabel)

    # resize the colorbar
    pos.y0 += 0.022
    pos.y1 -= 0.05
    # pos.x1 += 0.10
    # arrange and adjust the position of each axis, ticks, and ticklabels
    ax1.set_position(pos)
    ax2.set_position(pos)
    ax1.xaxis.set_ticks_position('top') # set the position of the first axis to right
    ax1.xaxis.set_label_position('top') # set the position of the fitst axis to right
    ax1.set_xlabel('VR')
    ax2.xaxis.set_ticks_position('bottom') # set the position of the second axis to right
    ax2.xaxis.set_label_position('bottom') # set the position of the second axis to right
    # ax2.spines['left'].set_position(('outward', 50)) # adjust the position of the second axis
    ax2.set_xlabel(r'$\overline{VR}$')


def compute_stressdrop(M0, fc, cs):
    wc = 2 * np.pi * fc # r0 is defined with angular frequency (see Mclaskey et al. 2014; Udias et al., 2014)
    r0 = 2.34*cs/wc # we assume Brune source model
    del_sigma = (7/16) * (M0/(r0**3)) 
    return del_sigma, r0

def compute_cornerfreq_madariaga(M0, R, cs, wavetype):
    """
    return corner frequency using Madariaga's crack model
    # note: we assume cp/cs = 1.732
    """
    if wavetype=="P":
        k = 0.32
    elif wavetype=="S":
        k = 0.21
    else:
        print("error: wavetype is either P or S.")
        return None

    return k * cs / R


def compute_scaling_M0(fc, delsig, cs):
    """
    return M0 to plot scaling line of stress drop
    """
    return (16/7) * delsig * (((2.34*cs)/(2*np.pi*fc))**3)


