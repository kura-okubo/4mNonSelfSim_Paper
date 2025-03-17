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
from matplotlib import ticker, cm
from matplotlib.colors import LogNorm
from matplotlib import ticker, cm, colors
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from compute_CF_func import *


def store_trace(datmat, fs, st_stats, ev_t, Nsensor = 32):

	# store data to trace
	st = Stream()
	for stid in tqdm(range(Nsensor)): 
		tr = Trace()
		tr.data = datmat[:, stid].view()
		tr.stats.network = "FB"
		tr.stats.station = "OL{:02d}".format(stid+1)
		tr.stats.channel = "Z"
		tr.stats.sampling_rate = fs
		tr.stats.starttime=st_stats + timedelta(seconds=float(ev_t))
		st.append(tr)

	return st

#----------------------------------------------------------#
def gui_pick_arrival(fig, st_obs):

	# decrare global variable
	global trace_id, tpick, tpick_trial, h_vl, ax1, iter_pick

	"""
	GUI to pick first arrival
	"""

	titlestr = "'a' pick 'c' compute location 'u' update arrival time 'w' save 'd' delete 'q' exit."

	def init_tpick():
		global tpick
		tpick = {}
		for i, tr in enumerate(st_obs):
			stnm = tr.stats.station
			# initialize tpick with large number; if arrival is not chosen, they are sorted out.
			tpick[stnm] = 1e2

	def find_nearest(array, value):
		array = np.asarray(array)
		idx = (np.abs(array - value)).argmin()
		return array[idx]

	def plot_trace(trace_id):
		global ax1, h_vl
		
		t = st_obs[0].times()*1e3
		plt.clf()
		tr = st_obs[trace_id]
		ax1 = fig.add_subplot(gs_1[:, :])

		ax1.plot(tr.times()*1e3, tr.data/np.max(np.abs(tr.data)), "k-") # we noramlize the amplitude for the sake of picking
		# ax1.plot(tr.times()*1e3, tr.data, "k-")

		# Comment out the autpick function
		# ax2 = ax1.twinx()
		# fap = st_obs[trace_id].fap 
		# ax2.plot(t, fap, "b:", alpha=0.5)
		# ax2.axhline(autopick_threshold, c="k", ls="--")

		ax1.set_xlabel('Time [ms]')
		ax1.set_ylabel('Amplitude', color='k')
		#   ax2.set_ylabel('Characteristic function', color='b')
		#   ax1.set_xlim([t[t_pick0]-0.1, t[t_pick0]+0.1])
		ax1.set_xlim([t[0], t[-1]])
		# ax1.set_ylim([-10.0, 10.0])
		# ax1.set_ylim([-5.0, 5.0])
		ax1.set_ylim([-1, 1])

		# check and get coordinate
		stnm = tr.stats.station
		ax1.set_title("{}\nevent{:04d} {} {} arrival time: {:4.5f} [ms]".format(titlestr, st_obs.gougeevent_id, st_obs.event_type, stnm, tpick[stnm]))
		h_vl = ax1.axvline(tpick[stnm], c="g", ls="-");
		fig.tight_layout()

	def plot_moveout(trace_id):
		t = st_obs[0].times()*1e3
		ax3 = fig.add_subplot(gs_2[:, :])
		ax3.clear()
		plotAEnorm = 2e-3
		plot_bar_length = 50
		skip_step = 1 #20


		for i, tr in enumerate(st_obs):
			stnm = tr.stats.station
			st_x = st_obs.channel_loc[stnm][0]
			if i == trace_id:
				lc = "b"
			else:
				lc = "k"

			ax3.plot(t[::skip_step], tr.data[::skip_step]/plotAEnorm + st_x, "-", c=lc) # decimate to fasten the plotting
			ax3.plot([tpick[stnm], tpick[stnm]], [st_x-plot_bar_length/2, st_x+plot_bar_length], "r-", zorder=100)

		ax3.set_xlabel('Time [ms]')
		ax3.set_ylabel('Distance [mm]', color='k')
		ax3.set_xlim([t[0], t[-1]])
		ax3.set_ylim([-50, 4150])
		fig.tight_layout()

	def on_press_all(event):
		global ax1, tpick, trace_id, iter_pick

		t = st_obs[0].times()*1e3
		stnm = st_obs[trace_id].stats.station

		if event.key == 'a': # modify pick time
			h_vl.set_color('g')			
			tpick[stnm] = find_nearest(t, event.xdata) # update pick time
			ax1.set_title("{}\nevent{:04d} {} {} arrival time: {:4.5f} [ms]".format(titlestr, st_obs.gougeevent_id, st_obs.event_type, stnm, tpick[stnm]))
			h_vl.set_xdata(tpick[stnm])
			plot_moveout(trace_id)
			plt.draw()

		elif event.key in ["up", "right"]: # move station id
			if trace_id < 16:
				trace_id = trace_id + 16
			elif trace_id <= 30:
				trace_id = trace_id - 15
			else: # loop trace id
				trace_id = 0
			plot_trace(trace_id)
			plot_moveout(trace_id)
			plt.draw()

		elif event.key in ["down", "left"]: # move station id
			if trace_id == 0:
				trace_id = 31
			elif trace_id < 16:
				trace_id = trace_id + 15
			else: # loop trace id
				trace_id = trace_id - 16
			plot_trace(trace_id)
			plot_moveout(trace_id)
			plt.draw()

		elif event.key == "c":
			global tpick_trial

			"""
			compute source location
			"""

			ax1.set_title("{}\nevent{:04d} {} {} arrival time: {:4.5f} [ms] location done!".format(titlestr, st_obs.gougeevent_id, st_obs.event_type, stnm, tpick[stnm]))

			datacase = st_obs.runID+"__{:04d}".format(st_obs.gougeevent_id)
			Vrng = [6200] # estimated P wave
			gridsize = 0.25e-3
			vfixed_id = 0
			PickNumStation = 4
			channel_loc = st_obs.channel_loc

			T, X, Y, V, R, Ermat, Vmat, Tmat, Xrng, Yrng  = eval_biax_4m_eventloc(tpick, st_obs, channel_loc, Vrng, PickNumStation=PickNumStation, 
                                                                  gridsize=gridsize, datacase=datacase, vfixed_id=vfixed_id)
			# compute synthetick p arrival
			tpick_trial = {}
			for stnm, _ in tpick.items():
				xloc = channel_loc[stnm][0]/1e3
				yloc = channel_loc[stnm][1]/1e3
				zloc = channel_loc[stnm][2]/1e3
				dist=np.sqrt((xloc-X)**2+(yloc-Y)**2+zloc**2)
				tpick_trial[stnm] = T+dist/V * 1e3 #[ms]

			## Compute distance with inverted location

			# AEsensor_newloc={}
			# source_loc=[X*1e3, Y*1e3, 0]

			# for i, stnm in enumerate(st_obs.channel_loc.keys()):
			#     station_loc = np.array(st_obs.channel_loc[stnm])
			#     dist = np.linalg.norm(source_loc-station_loc)
			#     AEsensor_newloc[stnm] = [dist, station_loc[0], station_loc[1], station_loc[2]]


			vmin = 5e-1 #5e-1
			vmax = 1e1 #2e1
			cmap = cm.jet

			Ermat = np.ma.masked_where(Ermat <= 0, Ermat)

			lev_exp = np.arange(np.floor(np.log10(vmin))-1e-9, np.ceil(np.log10(vmax))+1e-9, 0.01)
			levs = np.power(10, lev_exp)

			ax4 = fig.add_subplot(gs_3[:, :])
			ax4.clear()

			cs = ax4.contourf(Xrng*1e3, Yrng*1e3, np.transpose(1e3*np.sqrt(Ermat)), levs, 
			                 norm=LogNorm(vmin=vmin,vmax=vmax), cmap=cmap, extend='both') # plot in micro second

			ax4.set_xlabel("x [mm]")
			ax4.set_ylabel("y [mm]")
			ax4.invert_yaxis()
			ax4.set_aspect('equal', adjustable='box')

			cticks = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2]

			sm = plt.cm.ScalarMappable(cmap=cmap, norm=LogNorm(vmin=vmin,vmax=vmax))
			cbar = plt.colorbar(sm, ax=ax4, orientation='vertical', ticks=cticks, pad=0.05, extend='both')

			cbar.set_label(r'STD pick time diff. [$\mu$s]')

			# decorate figures
			# plot best location 
			ax4.plot(X*1e3, Y*1e3, marker='x', markersize=10, color='w', mew=1)

			# plot location of stations
			tpick_sorted = sorted(tpick.items(), key=lambda x:x[1], reverse=False)
			xlimlist = []
			for i, (stnm, t1) in enumerate(tpick_sorted[:PickNumStation]):
			    if np.isnan(t1):
			        # skip this station
			        continue;

			    x1 = channel_loc[stnm][0]
			    y1 = channel_loc[stnm][1]
			    z1 = channel_loc[stnm][2]
			    xlimlist.append(x1)
			    # upper side
			    if y1 >0:
			        ax4.plot(x1, y1-4, marker='^', markersize=10, color='c', mec='k')
			    elif y1 < 0:
			        ax4.plot(x1, y1+4, marker='v', markersize=10, color='c', mec='k')

			xmin = np.sort(xlimlist)[0]-10
			xmax = np.sort(xlimlist)[-1]+10
			ax4.set_xlim(xmin, xmax)
			ax4.set_yticks([-50, -25, 0, 25, 50])

			pos1 = ax4.get_position() # get the original position 
			pos2 = [pos1.x0 + 0.01, pos1.y0 + 0.05,  pos1.width, pos1.height*0.95] 
			ax4.set_position(pos2) # set a new position
			# fig.tight_layout()

			iter_pick = iter_pick+1
			figname = st_obs.eventoutdir+"/errvar_contour__{}__{:04d}__iter__{:02d}.png".format(st_obs.runID, st_obs.gougeevent_id, iter_pick)
			plt.savefig(figname, dpi=150)

			# keep temporal file to store location and origin time
			origin_abs_t = st_obs.init_abs_t + T*1e-3 # T is in ms.
			tpick_rel = {}
			for key in tpick:
				tpick_rel[key] = ( st_obs.init_abs_t + tpick[key]*1e-3 ) - origin_abs_t

			eventinfo = {
						"init_abs_t": st_obs.init_abs_t,
						"origin_abs_t": origin_abs_t, 
						"T": T,
						"X": X,
						"Y": Y,
						"V": V,
						"R": R,
						"tpick_rel":tpick_rel,
						}

			print(eventinfo["init_abs_t"], eventinfo["origin_abs_t"], eventinfo["T"])
			with open(fotmpeventname, 'wb') as fo:
				pickle.dump(eventinfo, fo, protocol=pickle.HIGHEST_PROTOCOL)
			ax1.set_title("{}\nevent{:04d} {} {} arrival time: {:4.5f} [ms] t_origin:{:4.2f}".format(titlestr, st_obs.gougeevent_id, st_obs.event_type, stnm, tpick[stnm], T))

		elif event.key == "u":
			"""
			update tpick with synthetic arrival time
			"""
			tpick = tpick_trial
			plot_trace(trace_id)
			plot_moveout(trace_id)
			plt.draw()

		elif event.key == "w":
			"""
			save tpick into file
			"""
			with open(fopickname, 'wb') as fo:
				pickle.dump(tpick, fo, protocol=pickle.HIGHEST_PROTOCOL)
			ax1.set_title("{}\nevent{:04d} {} {} arrival time: {:4.5f} [ms] save tpick done!".format(titlestr, st_obs.gougeevent_id, st_obs.event_type, stnm, tpick[stnm]))
			plt.draw()
			# rename event loc file
			os.rename(fotmpeventname, foeventname)

		elif event.key == "d":
			"""
			discard this event
			"""
			if os.path.exists(fopickname):
				os.remove(fopickname)
			if os.path.exists(fotmpeventname):
				os.remove(fotmpeventname)
			if os.path.exists(foeventname):
				os.remove(foeventname)

			ax1.set_title("{}\nevent{:04d} {} {} arrival time: {:4.5f} [ms] pickfile removed.".format(titlestr, st_obs.gougeevent_id, st_obs.event_type, stnm, tpick[stnm]))
			plt.draw()

		elif event.key == "x":
			"""
			exit the GUI
			"""
			print("exit GUI with x key.")
			exit()

	#--------#
	gs_master = GridSpec(nrows=3, ncols=1, height_ratios=[1, 1.5, 1.5])
	gs_1 = GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=gs_master[0, 0])
	gs_2 = GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=gs_master[1, 0])
	gs_3 = GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=gs_master[2, 0])

	# output name of tpick
	fopickname = os.path.join(st_obs.eventoutdir, "tpick__{}__{:04d}.pickle".format(st_obs.runID, st_obs.gougeevent_id))
	fotmpeventname = os.path.join(st_obs.eventoutdir, "tmp_eventloc__{}__{:04d}.pickle".format(st_obs.runID, st_obs.gougeevent_id))
	foeventname = os.path.join(st_obs.eventoutdir, "eventloc__{}__{:04d}.pickle".format(st_obs.runID, st_obs.gougeevent_id))

	# initialize canvas with event
	cid = fig.canvas.mpl_connect('key_press_event', on_press_all) # turn on interactive events

	# initialize tpick arrival time
	if os.path.exists(fopickname):
		print("tpick file found: {}".format(fopickname))
		with open(fopickname, 'rb') as fi:
			tpick = pickle.load(fi)

	else:
		init_tpick()

	# plot figures
	trace_id = 0
	iter_pick = 0
	plot_trace(trace_id)
	plot_moveout(trace_id)
	plt.show()

