import os
import obspy
from obspy import read, Stream, Trace
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime
from datetime import timedelta
from tqdm import tqdm
import warnings
import matplotlib as mpl
from scipy.interpolate import LSQUnivariateSpline

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 12
os.environ['TZ'] = 'GMT' # change time zone to avoid confusion in unix_tvec conversion

import xml.etree.ElementTree as ET

def M02Mw(M0):
    return (np.log10(M0) - 9.1) * 2.0 / 3.0 # synchronized with OpenSWPC : moment_magnitude ( m0 )

def Mw2M0(Mw):
    return 10**( 1.5 * Mw + 9.05) # synchronized with OpenSWPC : seismic_moment ( mw )
   
def read_hsu1985_crossverification(finame):
    """
    parse xml struct obtained from MATLAB software tool of generalized ray theory
    """
    tree = ET.parse(finame)
    root = tree.getroot()

    npts = int(root.find("npoint").text)
    tvec = np.array([float(x.text) for x in root.findall("tvec")])
    dt = float(root.find("dt").text)

    st = Stream()
    # roop for stations
    for stationname in ["R01", "R02", "R03"]:
        # roop for components
        for gcompname in ["G11", "G22", "G33","G111","G112","G113","G121","G122","G123","G131","G132","G133",
                          "G211","G212","G213","G221","G222","G223","G231","G232","G233",  
                            "G311","G312","G313","G321","G322","G323","G331","G332","G333"]:

            
            comp = root.find(stationname).find(gcompname)
            if comp is None:
                # print("{} not exist. skipping.".format(gcompname))
                continue

            tr = Trace()
            tr.stats.network = "HSU"
            tr.stats.station = stationname
            tr.stats.channel = gcompname
            tr.stats.starttime = datetime.datetime(2022,10,18,8)
            tr.stats.delta = dt
            tr.data = np.array([float(x.text) for x in root.find(stationname).find(gcompname).findall("data")])
            st.append(tr)
            
    return st


def stf_cosine(t, TR, fz):
    '''
        source time function of cosine wavelet
        https://tktmyd.github.io/OpenSWPC/English/2._Parameter_Settings/0207_source/
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




def stf_herzian_mclaskey2009(t, rho, R, v, E1, nu1, E2, nu2):
    '''
        source time function of Herzian solution
    '''
    def get_delta(E1, mu1, E2, mu2):
        del1 = (1-nu1**2)/(np.pi*E1)
        del2 = (1-nu2**2)/(np.pi*E2)
        return (del1, del2)
    
    
    def get_contact_time(rho, del1, del2, R, v):
        return 4.53*(( 4*rho  * np.pi * (del1 + del2) /3 )**(2/5)) * R * (v **(-1/5))

    def get_maximum_force(rho, del1, del2, R, v):
        return 1.917 * (rho**(3/5)) * ((del1 + del2)**(-2/5)) * (R**2) * (v**(6/5))

        
    del1, del2 = get_delta(E1, nu1, E2, nu2)
    
    # print(del1, del2, rho, R, v)
    tc = get_contact_time(rho, del1, del2, R, v)
    fmax = get_maximum_force(rho, del1, del2, R, v)
    print(tc, fmax)
    
    stf = np.zeros(len(t))
    for i, tt in enumerate(t):
        if 0<tt and tt < tc:
            stf[i] = fmax*(np.sin(np.pi*tt/tc))**(3/2)
        else:
            stf[i] = 0
            
    return stf




def compute_synthetic_waveform_vz_openswpc(M0, mij, st_syn_sta):
    """
    compute synthesized waveform
    """
    if len(st_syn_sta) != 6:
        warnings.warn("number of green's tensor is not 6.") 
    Gzxx = st_syn_sta.select(channel="G_Vz_mxx")[0].data
    Gzyy = st_syn_sta.select(channel="G_Vz_myy")[0].data
    Gzzz = st_syn_sta.select(channel="G_Vz_mzz")[0].data
    Gzxy = st_syn_sta.select(channel="G_Vz_mxy")[0].data
    Gzxz = st_syn_sta.select(channel="G_Vz_mxz")[0].data
    Gzyz = st_syn_sta.select(channel="G_Vz_myz")[0].data
    
    vz = M0 * (mij[0, 0] * Gzxx + mij[1, 1] * Gzyy + mij[2, 2] * Gzzz +
               mij[0, 1] * Gzxy + mij[0, 2] * Gzxz + mij[1, 2] * Gzyz)
    
    return vz

def compute_synthetic_waveform_vx_hsu1985(mij, stf, st_hsu):
    """
    compute synthesized waveform associated with analytical solution by Hsu1985
    """

    if len(st_hsu) != 9:
        warnings.warn("number of green's tensor is not 9.") 
        
    ## Assert the symmetricity of Green's tensor
    N = st_hsu.select(channel="G111")[0].stats.npts
    dt = st_hsu.select(channel="G111")[0].stats.delta
    vx_sum = np.zeros((N,))
    
    for i in range(3):
        for j in range(3):
            if "{}{}".format(i+1, j+1) in ["12", "21", "23", "32"]:
                continue; 

            green_tmp = st_hsu.select(channel="G1{}{}".format(i+1, j+1))[0].data
#             print(st_hsu.select(channel="G3{}{}".format(i+1, j+1))[0].stats)
            # convolve the source time function
            vx_sum += mij[i, j] * np.convolve(green_tmp, stf, mode='same') * dt # Note that M0 is already multiplied in the stf.          

            
    return vx_sum

def compute_synthetic_waveform_vy_hsu1985(mij, stf, st_hsu):
    """
    compute synthesized waveform associated with analytical solution by Hsu1985
    """

    if len(st_hsu) != 9:
        warnings.warn("number of green's tensor is not 9.") 
        
    ## Assert the symmetricity of Green's tensor
    N = st_hsu.select(channel="G211")[0].stats.npts
    dt = st_hsu.select(channel="G211")[0].stats.delta
    vy_sum = np.zeros((N,))
    
    for i in range(3):
        for j in range(3):
            if "{}{}".format(i+1, j+1) in ["11", "13", "22", "31", "33"]:
                continue;

            green_tmp = st_hsu.select(channel="G2{}{}".format(i+1, j+1))[0].data
#             print(st_hsu.select(channel="G3{}{}".format(i+1, j+1))[0].stats)
            # convolve the source time function
            vy_sum += mij[i, j] * np.convolve(green_tmp, stf, mode='same') * dt # Note that M0 is already multiplied in the stf.          

            
    return vy_sum


def compute_synthetic_waveform_vz_hsu1985(mij, stf, st_hsu):
    """
    compute synthesized waveform associated with analytical solution by Hsu1985
    """

    if len(st_hsu) != 9:
        warnings.warn("number of green's tensor is not 9.") 
        
    ## Assert the symmetricity of Green's tensor
    N = st_hsu.select(channel="G311")[0].stats.npts
    dt = st_hsu.select(channel="G311")[0].stats.delta
    vz_sum = np.zeros((N,))
    
    for i in range(3):
        for j in range(3):
            if "{}{}".format(i+1, j+1) in ["12", "21", "23", "32"]:
                continue;

            green_tmp = st_hsu.select(channel="G3{}{}".format(i+1, j+1))[0].data
#             print(st_hsu.select(channel="G3{}{}".format(i+1, j+1))[0].stats)
            # convolve the source time function
            vz_sum += mij[i, j] * np.convolve(green_tmp, stf, mode='same') * dt # Note that M0 is already multiplied in the stf.          

            
    return vz_sum


def compute_synthetic_waveform_vz_herrmann96(phi_deg, M0, mij, stf, st_syn_herrmann):
    """
    compute synthesized waveform associated with Computer program in Seismology
    Herrmann (2013)
    Note: As the formulation of uz (vz) explicitly includes the azimuth (see also Jost and Hermann, 1989)
    so the azimuth could be balanced with the moment tensor solution.
    """

    phi = np.deg2rad(phi_deg)

    if len(st_syn_herrmann) != 15:
        warnings.warn("number of green's tensor with Herrmann96 is not 15.") 

    ZDD = st_syn_herrmann.select(channel="ZDD")[0].data
    ZSS = st_syn_herrmann.select(channel="ZSS")[0].data
    ZEX = st_syn_herrmann.select(channel="ZEX")[0].data
    ZDS = st_syn_herrmann.select(channel="ZDS")[0].data

    # 1. Synthesize the time history of Green's function
    # NOTE: multiply 2 for Mxy and Myx
    
    green_tmp = mij[0, 0] * ((ZSS/2) * np.cos(2*phi) - ZDD/6 + ZEX/3) \
                + mij[1, 1] * (-(ZSS/2) * np.cos(2*phi) - ZDD/6 + ZEX/3) \
                + mij[2, 2] * (ZDD/3 + ZEX/3) \
                + 2*mij[0, 1] * (ZSS * np.sin(2*phi)) \
                + 2*mij[0, 2] * (ZDS * np.cos(phi)) \
                + 2*mij[1, 2] * (ZDS * np.sin(phi))
    

    N = st_syn_herrmann.select(channel="ZDD")[0].stats.npts
    dt = st_syn_herrmann.select(channel="ZDD")[0].stats.delta

    # convolve source time function
    vz = M0*np.convolve(green_tmp, stf, mode='same') * dt

    return vz



def compute_momenttensor_lab(strike_deg, rake_deg):
    '''
    Return moment tensor associated with events on 4m biax fault
    assuming that dip = 90, strike = 0 in the simulation coordinates (sidecoord: x along fault, y along side of fault, z upwards.) 
    input:
        rake (radian): rake of slip
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
    mij[1][2] = np.sin(rake)*np.cos(strike) # positive in southern side of specimen; i.e. positive towards upward zenith
    mij[2][0] = mij[0][2]
    mij[2][1] = mij[1][2]
    
    mag = np.linalg.norm(mij, ord=None)
    
    return mij/mag



