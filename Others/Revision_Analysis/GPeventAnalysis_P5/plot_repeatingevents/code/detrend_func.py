import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy import integrate, signal


def STF_multiwindow_detrend(tvec, data_raw, onset, polyfit_winlens, polyfit_Rcontinu_buffer, 
                            smooth_lowpass_freqmax, polyords=[3,2,3], debugplot=False, eventdata=None, max_pwinlen=5e-6):
    
    """
    Apply the detrend by the polynomial fitting on the backward, STF and the P coda time windows
    to remove the long-period noise from the source time function.

    The data is separted by the three time windows as follows:
    
    | backward | STF | P coda |
    variable: onset  STF_Rcontinu
            
    Process flow:
    1. Compute polyfit on the backward and the STF time windows.
    2. Merge the fit functions as the trend trace at the continuity between backward and the STF specified by the `onset`.
    3. Find the best continuity point of the trend trace between STF and Pcoda, where the residuals
    between the merged trend and the data is small.
    4. Compute polyfit on the P coda window.
    5. Merge the fitting function.
    6. Apply the lowpass filter to smooth the discontinuity at the continuities.
    7. Remove the trend from the data.

    Arguments:

    tvec : [float vector] time vector
    data_raw : [float vector] waveform data
    onset :[float] onset time of the P wave pulse
    polyfit_winlens: [polyfitwinlen_backward, polyfitwinlen_STF, polyfitwinlen_Pcoda][float] the boundaries to
    apply the polynomial fitting.
    polyfit_Rcontinu_buffer: [float] search the STF_Rcontinu from the polyfitwinlen_STF - STF_Rcontinu
    smooth_winlen: [float] window length of the moving window smoothing
    polyords: [int vector] the orders of the polynomial fitting.
    """

    data = copy.deepcopy(data_raw) # copy the raw data to avoid unexpected change
    Ndata = len(data)
    dt = tvec[1]-tvec[0]

    onset_k = int(onset/dt)
    
    polyfitwinlen_backward, polyfitwinlen_STF, polyfitwinlen_Pcoda = polyfit_winlens
        
    polyfit_backward_k = int(polyfitwinlen_backward/dt)
    polyfit_STF_k = int(polyfitwinlen_STF/dt)
    polyfitwinlen_Pcoda_k = int(polyfitwinlen_Pcoda/dt)

    w_backward = np.zeros(Ndata)
    w_STF = np.zeros(Ndata)
    w_Pcoda = np.zeros(Ndata)
    tr_detrend_merged = np.zeros(Ndata)

    # 1. Polyfit on the backward and the STF time windows
    w_backward[onset_k-polyfit_backward_k:onset_k] = 1
    w_STF[onset_k:onset_k+polyfit_STF_k] = 1
    
    # weight at the continuity
    w_backward[onset_k] = 1e3
    w_STF[onset_k] = 1e3

    k_poly_backward = np.polynomial.polynomial.polyfit(tvec, data, polyords[0], w=w_backward)
    k_poly_STF = np.polynomial.polynomial.polyfit(tvec, data, polyords[1], w=w_STF)

    tr_poly_backward = np.polynomial.polynomial.polyval(tvec, k_poly_backward)
    tr_poly_STF = np.polynomial.polynomial.polyval(tvec, k_poly_STF)

    # 2. Merge the fit function
    tr_detrend_merged[:onset_k] = tr_poly_backward[:onset_k]
    tr_detrend_merged[onset_k:] = tr_poly_STF[onset_k:]

    # 3. Find the continuity between the data and the merged trend trace
    residu = data - tr_detrend_merged
    res_eps = np.abs(residu[onset_k:onset_k+polyfit_STF_k]).max() * 1e-3 # find the points where the residual is less than this eps.
    cross_all = np.where(np.abs(residu) < res_eps)[0]

    # pick the time closest to the boundary between the STF and Pcoda
    cross_STFwin_closest_k = np.argmin(np.abs(tvec[cross_all] - (onset+polyfitwinlen_STF)))
    STF_Rcontinu_k = cross_all[cross_STFwin_closest_k]
    STF_Rcontinu_t = tvec[STF_Rcontinu_k]
    tdiff_Rcontinu = np.abs(STF_Rcontinu_t - (onset+polyfitwinlen_STF))

    # threshold the distance from the (onset+polyfitwinlen_STF). If the difference is too large,
    # we set the Rcontinu at the polyfit_Rcontinu_buffer as it is not expected to set the Rcontinu
    # during the P wave pulse
    if tdiff_Rcontinu > polyfit_Rcontinu_buffer:
        print("tdiff of Rcontinu is larger than polyfit_Rcontinu_buffer: set the Rcontinu at the polyfit_Rcontinu_buffer.")
        STF_Rcontinu_k = np.where(tvec>(onset+polyfitwinlen_STF-polyfit_Rcontinu_buffer))[0][0]
        STF_Rcontinu_t = tvec[STF_Rcontinu_k]
    
    # 4. Compute polyfit on the P coda window.
    polyfit_Pcoda_k = int((polyfitwinlen_Pcoda)/dt)
    w_Pcoda[STF_Rcontinu_k:STF_Rcontinu_k+polyfit_Pcoda_k] = 1
    # weight at the continuity
    w_Pcoda[STF_Rcontinu_k] = 1e3
    
    k_poly_Pcoda = np.polynomial.polynomial.polyfit(tvec, data, polyords[2], w=w_Pcoda)

    tr_poly_Pcoda = np.polynomial.polynomial.polyval(tvec, k_poly_Pcoda)

    # 5. Merge the fitting function
    tr_detrend_merged[STF_Rcontinu_k:] = tr_poly_Pcoda[STF_Rcontinu_k:]
    
    # 6. Apply the two-way Low-pass filter (we deprecated moving window smoothing)
    # smooth_winlen_k = int(smooth_winlen/dt)
    # #smooth after computing Pcoda fit
    # tr_detrend_merged_smoothed = np.convolve(tr_detrend_merged, np.ones(smooth_winlen_k)/smooth_winlen_k, mode='same')

    b_trendsmooth, a_trendsmooth = signal.butter(3, smooth_lowpass_freqmax, 'lowpass', fs=(1/dt), output='ba')
    tr_detrend_merged_smoothed = signal.filtfilt(b_trendsmooth, a_trendsmooth, tr_detrend_merged, method='gust')
    
    # 7. Remove the the trend
    data_detrend_multiwindows = data - tr_detrend_merged_smoothed

    # if removeoffset:
    #     # remove the linear trend, and the constant offset on the backward window
    #     # 1. remove the linear trend on the backward window
    #     k_poly_backward_lintrend = np.polynomial.polynomial.polyfit(tvec, data_detrend_multiwindows, 1, w=w_backward)
    #     tr_poly_backward_lintrend = np.polynomial.polynomial.polyval(tvec, k_poly_backward_lintrend)
    #     print(tr_poly_backward_lintrend[onset_k], data_detrend_multiwindows[onset_k])
    #     # make continuity at onset
    #     tr_detrend_backward_lintrend = np.ones(Ndata) * tr_poly_backward_lintrend[onset_k]
    #     tr_detrend_backward_lintrend[:onset_k] = tr_poly_backward_lintrend[:onset_k]
    #     # tr_detrend_backward_lintrend_smoothed =  np.convolve(tr_detrend_backward_lintrend, np.ones(smooth_winlen_k)/smooth_winlen_k, mode='same')
    #     tr_detrend_backward_lintrend_smoothed = tr_detrend_backward_lintrend
    #     data_detrend_multiwindows[:onset_k] -= tr_detrend_backward_lintrend_smoothed[:onset_k] - tr_poly_backward_lintrend[onset_k]

    #     # 1. remove the constant offset
    #     st_offset = int((onset-polyfit_winlens[0])/dt)
    #     et_offset = int(onset/dt)
    #     constant_offset = np.mean(data_detrend_multiwindows[st_offset:et_offset])
    #     data_detrend_multiwindows -= constant_offset
        
    # else:
    #     # set for the debug plot
    #     tr_detrend_backward_lintrend_smoothed = np.zeros(Ndata)
    #     constant_offset = 0

    # tr_detrend_backward_lintrend_smoothed = np.zeros(Ndata)
    # constant_offset = 0
    
    if debugplot:
        debugplot_STF_multiwindow_detrend(tvec, data_raw, data_detrend_multiwindows,
                                          w_backward, w_STF, w_Pcoda,
                                          tr_poly_backward, tr_poly_STF, tr_poly_Pcoda,
                                          tr_detrend_merged, tr_detrend_merged_smoothed,
                                          #tr_detrend_backward_lintrend_smoothed, constant_offset,
                                          onset, STF_Rcontinu_k, STF_Rcontinu_t, eventdata, max_pwinlen)
        
    return data_detrend_multiwindows




def debugplot_STF_multiwindow_detrend(tvec_raw, data_raw,  data_detrend_multiwindows,
                                          w_backward, w_STF, w_Pcoda,
                                          tr_poly_backward, tr_poly_STF, tr_poly_Pcoda,
                                          tr_detrend_merged, tr_detrend_merged_smoothed,
                                          # tr_detrend_backward_lintrend_smoothed, constant_offset,
                                          onset, STF_Rcontinu_k, STF_Rcontinu_t, eventdata, max_pwinlen=5e-6):

    tvec = copy.deepcopy(tvec_raw)*1e6
    dt = tvec[1]-tvec[0]
    onset_k = int(onset/(dt*1e-6))

    fig, axs = plt.subplots(3, 1, figsize=(8,7))

    # 1. Raw trace and the polynomial fitting functions
    # scaling by eventdata["k_M0uz_event"] 
    axs[0].plot(tvec, data_raw/eventdata["k_M0uz_event"]*1e3, "k-", label="before detrend", zorder=3)
    axs[0].plot(tvec, tr_poly_backward/eventdata["k_M0uz_event"]*1e3, "b:", label=None, zorder=-1)
    axs[0].plot(tvec, tr_poly_STF/eventdata["k_M0uz_event"]*1e3, "g:", label=None, zorder=-1)
    axs[0].plot(tvec, tr_poly_Pcoda/eventdata["k_M0uz_event"]*1e3, "y:", label=None, zorder=-1)
    # axs[0].plot(tvec, tr_detrend_merged/eventdata["k_M0uz_event"]*1e3, "--", c="crimson", label="trend trace", zorder=4)
    axs[0].plot(tvec, tr_detrend_merged_smoothed/eventdata["k_M0uz_event"]*1e3, "-", c="r", label="merged trend", zorder=4)

    # plot window boundaries
    axs[0].axvline(onset*1e6, c="gray", ls="--", lw=0.75)
    axs[0].axvline(STF_Rcontinu_t*1e6, c="gray", ls="--", lw=0.75)
    axs[0].plot(tvec[STF_Rcontinu_k], tr_detrend_merged_smoothed[STF_Rcontinu_k]/eventdata["k_M0uz_event"]*1e3, "o",
                zorder=5, mfc="w", mec="k", ms=4)
    
    ylimit_vel1 = 1.2*np.max(np.abs(data_raw))*np.array([-1, 1])/eventdata["k_M0uz_event"]*1e3
    axs[0].set_ylim(ylimit_vel1)

    axs[0].set_ylabel("Velocity [mm/s]")
    axs[0].legend(loc=2)

    # plot the polyfit weigh 
    ax2 = axs[0].twinx()
    ax2.plot(tvec, w_backward, c="orange", ls=":", lw=0.75)
    ax2.plot(tvec, w_STF, c="cyan",ls=":", lw=0.75)
    ax2.plot(tvec, w_Pcoda, c="magenta", ls=":", lw=0.75)
    ax2.set_ylim([-1, 2])
    ax2.set_ylabel("Polyfit weight")
    
    # Comparison before and after the detrend
    axs[1].plot(tvec, data_raw/eventdata["k_M0uz_event"]*1e3, label="before detrend", c="k", zorder=1);
    axs[1].plot(tvec, data_detrend_multiwindows/eventdata["k_M0uz_event"]*1e3, c="r", label="after detrend", zorder=2);
    
    # # plot additional detrend on the backward window
    # axs[1].plot(tvec, tr_detrend_backward_lintrend_smoothed/eventdata["k_M0uz_event"]*1e3, label=None, c="b", ls=":", lw=0.75, zorder=1);
    # axs[1].plot(0.4, constant_offset/eventdata["k_M0uz_event"]*1e3, label=None, c="g", marker="<", ms=6, zorder=1);
        
    axs[1].set_ylabel("Velocity [mm/s]")
    axs[1].legend(loc=2)

    ylimit_vel2 = 1.2*np.max(np.abs(data_detrend_multiwindows[onset_k:]))*np.array([-1, 1])/eventdata["k_M0uz_event"]*1e3
    axs[1].set_ylim(ylimit_vel2)

    # compute integral
    data_M0dot = integrate.cumulative_trapezoid(data_detrend_multiwindows, dx=dt*1e-6, initial=0)
    offset = data_M0dot[onset_k]
    data_M0dot -= offset

    ylimit_M0dot = 2.0*np.max(np.abs(data_M0dot[onset_k:onset_k+int(max_pwinlen/(dt*1e-6))]))*np.array([-1, 1])/1e6
    axs[2].plot(tvec, data_M0dot/1e6, c="r");
    axs[2].set_ylabel(r"$\dot{M}_0$ [MNm/s]")
    axs[2].set_ylim(ylimit_M0dot)
 
    # annotate tp and ts
    P_onset = eventdata["P_onset"]
    tp = eventdata["source_dist"]/eventdata["vp"] # [μs]
    ts = eventdata["source_dist"]/eventdata["vs"] # [μs]
    del_tpts = ts-tp # [μs]
  
    ty = axs[2].get_ylim()[1] * 0.7
    axs[2].plot(P_onset*1e6, ty, "kv")
    axs[2].text(P_onset*1e6, ty*1.1, "P", ha="center", clip_on=True)
    axs[2].plot((P_onset+del_tpts)*1e6, ty, "ks")
    axs[2].text((P_onset+del_tpts)*1e6, ty*1.1, "S", ha="center", clip_on=True)
    
    for i in range(3):
        axs[i].set_xlim([tvec[0], tvec[-1]])
        axs[i].set_xlabel("Time [μs]")

    axs[0].set_title(f'{eventdata["datacase"]} {eventdata["stnm"]}')
    plt.tight_layout()

    if "figname" in eventdata.keys():
        plt.savefig(eventdata["figname"], format="png", dpi=80)
        
    return fig