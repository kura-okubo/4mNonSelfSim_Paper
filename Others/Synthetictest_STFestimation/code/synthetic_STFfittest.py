import numpy as np
from scipy.interpolate import CubicSpline


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


def compute_res(x, STF_sensor, dt,  pwin_pre, residu_win=[0.5, 0.5], stf_type="cosine", debug=0):

    M0_try, TR_try, tshift_try = x
    
    """
    residu_win = [forward fraction, backward fraction]
    evaluate the residual from the peak to peak-residu_win[0]*TR, and to peak+residu_win[1]
    """
    #1. compute stf
    tvec_syn = np.linspace(0, TR_try, int(TR_try/dt))

    if stf_type=="kupper":
        STF_syn = stf_kupper(tvec_syn, TR_try, M0_try)
    elif stf_type=="cosine":
        STF_syn = stf_cosine(tvec_syn, TR_try, M0_try)
    else:
        raise ValueError("stf_type not found")

    # trim the synthetic STF
    st_syn = (0.5 - residu_win[0])*TR_try
    st_syn_k = int(st_syn/dt)
    et_syn_k = int((0.5 + residu_win[1])*TR_try/dt)
    N_res = et_syn_k-st_syn_k
    
    # trim the data STF
    st_obs_k = int((pwin_pre + tshift_try + (0.5-residu_win[0])*TR_try)/dt)
    et_obs_k = st_obs_k + N_res
    
    # onset_obs_k = int((pwin_pre + tshift_try)/dt)
    # assert onset_obs_k + N_res < len(STF_sensor)
    
    #4. compute residual
    if et_obs_k < len(STF_sensor):
        rmse = np.sqrt(np.mean((STF_sensor[st_obs_k:et_obs_k] - STF_syn[st_syn_k:et_syn_k])**2))
    elif st_obs_k > len(STF_sensor):
        print("debug data overflow")
        rmse = np.inf
        
    else:
        # print("debug data edge")
        et_obs_k = len(STF_sensor)
        Ndata = et_obs_k-st_obs_k
        # print(st_obs_k, et_obs_k,  len(STF_sensor), Ndata)
        rmse = np.sqrt(np.mean((STF_sensor[st_obs_k:et_obs_k] - STF_syn[st_syn_k:st_syn_k+Ndata])**2))
    
    if debug:
        return rmse, tvec_syn[st_syn_k:et_syn_k], STF_sensor[st_obs_k:et_obs_k], STF_syn[st_syn_k:et_syn_k]
    else:
        return rmse


def get_Qinv(freq, fq, Qinv):
    # interpolate the Q from the Qinv data
    cs = CubicSpline(fq, Qinv)
    
    Qinv_interp = np.zeros(len(freq))
    for i, ff in enumerate(freq):
        if ff<fq[0]:
            Qinv_interp[i] = Qinv[0] # extrapolate the minimum frequency Qinv
        elif ff>fq[-1]:
            Qinv_interp[i] = Qinv[-1] # extrapolate the maximum frequency Qinv
        else:
            Qinv_interp[i] = cs(ff)  
                
    return Qinv_interp