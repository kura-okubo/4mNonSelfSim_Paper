#!/usr/bin/env python
# coding: utf-8

# # AE locate events
# 
# 
# Locate the events with GUI picking.
# 
# - 2023.10.21 Kurama Okubo
# - 2024.10.30 update for the merged catalog.

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
import scipy.io as sio
import warnings
from obspy.core.utcdatetime import UTCDateTime    

from AE_compute_location_func import *

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 12
os.environ['TZ'] = 'GMT' # change time zone to avoid confusion in unix_tvec conversion
UTCDateTime.DEFAULT_PRECISION = 8 # increase the time precision

figdir = "../figure/AElocation"
if not os.path.exists(figdir):
    os.makedirs(figdir)

# Data outputdir
outdir = "../data/AElocation"
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Paramters
fs      = 10e6 #[Hz] Sampling frequency of data
Nsensor = 32 # Number of sensors
st_stats = UTCDateTime(2023, 5, 29) # starttime for the stats of trace: used to synchronize the data. Should be same with picking script.
expr_id= 87 # casename of experiment

ev_twinlen = 0.3e-3 #[s] window length of event
ev_pretrigger = 0.05e-3 #[s] duration of pretrigger

ev_twinlen_k = np.round(ev_twinlen*fs)

# Load AE location
channel_finame = '../../Others/AEchanneltable/AEsensorlocation_onFB03_table.csv'

df_array = pd.read_csv(channel_finame)

channel_loc={}

for i in range(len(df_array)):
    stnm = df_array.iloc[i].Instrument_Label
    xtemp = df_array.iloc[i].North.astype('float')
    ytemp = df_array.iloc[i].East.astype('float')
    ztemp = df_array.iloc[i].Down.astype('float')
    channel_loc[stnm] = [xtemp, ytemp, ztemp]
    
Nsensor = len(channel_loc)
# channel_loc

# Read picked catalog
# Read the csv pick datasheet
fi_catalog="../../Experiments/DetectEvent/data/p06_visual_pick_gougeevents_merged.csv"
columns = ["expr_id", "event_id", "event_loc", "picktime", "event_type", "rupturetype", "gougeevent_id", "doublecheck", "old_gougeevent_id"]
df_catalog = pd.read_csv(fi_catalog, skiprows=5, names=columns)
df_expr = df_catalog[(df_catalog["expr_id"]==f"fb03-{expr_id:03d}")]
df_expr.head()

# foreshock_eventset1 = [3, 16, 17, 30, 51, 52, 58, 61, 77, 91, 99, 106] # Ordinary, x=1750, rupturetype=1
# foreshock_eventset1 = [3, 16, 17, 30, 51, 52, 58, 61, 77, 91, 99, 106, 21, 32, 35, 67, 90, 97, 104] # Ordinary, x=1750, rupturetype=1 or 2
# foreshock_eventset2 = [3, 16, 17, 30, 51, 52, 58, 61, 77, 91, 99, 106, 21, 32, 35, 67, 90, 97, 104, \
            # 24, 36, 41, 43, 46, 70, 71, 81 ,82] # Ordinary, x=1750, rupturetype all (1, 2, and 3) # excluded 37 due to low S/N 
foreshock_eventset_merged = [4,   9,  18,  19,  20,  21,  24,  27,  30,  31,  37,  38,  40,
        43,  44,  49,  50,  52,  55,  59,  61,  62,  69,  72,  75,
        76,  77,  81,  85,  88,  89,  95,  99, 100, 102, 109, 110, 111,
       118, 120, 126, 128, 129, 131] # excluded 45 (old id:37) due to low S/N 

# foreshock_eventset_merged = [59,  61,  62,  69,  72,  75,
#         76,  77,  81,  85,  88,  89,  95,  99, 100, 102, 109, 110, 111,
#        118, 120, 126, 128, 129, 131]

for i, df_evt in df_expr.iterrows():
    # if not df_evt['gougeevent_id']==20:
    if not df_evt['gougeevent_id'] in foreshock_eventset_merged:
        continue; # skip events
    # df_evt = df_expr.iloc[12]
    df_evt

    # Load AE data from event mat data
    data_rootdir = f"/Volumes/4mGouge_WorkHDD/FB03data/4mBIAX_paper_tmp/p03_eventdata_FB03_{expr_id:03d}/"
    fname = f"eventdata_FB03_{expr_id:03d}_event{df_evt['event_id']:02d}"

    D = sio.loadmat(data_rootdir+fname)

    print(f"start processing eventid={df_evt['event_id']}: gougeevent_id={df_evt['gougeevent_id']:04d}")

    # read data
    ev_t = (df_evt["picktime"] - D["Tstart"])[0][0]
    init_abs_t = ev_t-ev_pretrigger #[s] absolute initial time of sta/lta trigger with pretrigger window
    print(ev_t, init_abs_t)
    init_k = int(init_abs_t*fs) # trim start time from picktime - pretrig
    datmat = D['AEdatmat'][init_k:int(init_k+ev_twinlen_k), :]

    # store data to stream
    st_obs = store_trace(datmat, fs, st_stats + df_evt["picktime"] - ev_pretrigger, ev_t)
    st_obs.channel_loc = channel_loc

    # store meta data
    outputdir = os.path.join(outdir, "arrivalpick")
    eventoutdir = outputdir+"/{:04d}".format(df_evt['gougeevent_id']) #i)
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    if not os.path.exists(eventoutdir):
        os.makedirs(eventoutdir)

    fig_snapshotdir = "../figure/AElocation/snapshot"
    if not os.path.exists(fig_snapshotdir):
        os.makedirs(fig_snapshotdir)

    st_obs.eventoutdir = eventoutdir
    st_obs.runID = f"fb03-{expr_id:03d}"
    st_obs.gougeevent_id = df_evt['gougeevent_id'] # index
    st_obs.event_type = df_evt["event_type"]
    st_obs.init_abs_t = D["Tstart"] + init_abs_t
    st_obs.fig_snapshotdir = fig_snapshotdir

    fig = plt.figure(figsize=(12, 10)) # initialize interactive figure

    # Apply bandpass filter
    # to detect Foreshocks with fc ~ 100kHz
    freq_min = 0.06e6 # [Hz]
    freq_max = 1.0e6 #0.6e6 #[Hz]
    st_obs.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=2, zerophase=True)

    # to detect low frequency event with fc ~ 20kHz
    # st_obs.filter('lowpass', freq=freq_max, corners=2, zerophase=True)

    gui_pick_arrival(fig, st_obs)

print(f"All events on fb03-{expr_id:03d} has been done!")

