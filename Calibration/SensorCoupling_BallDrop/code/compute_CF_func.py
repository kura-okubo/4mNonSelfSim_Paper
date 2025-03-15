import os
import numpy as np
from tqdm import tqdm
import warnings

def get_autopick_function(tr, FW, BW, st=0, et=np.inf):
    """
    return autopick function.
    
    index outside of BW and FW are zero-padded so that the data length of f is identical to original trace.
    """
    N = tr.stats.npts
    d = np.power(tr.data, 2)
    f = np.zeros(N)
    
    t = np.arange(0, (tr.stats.npts-1)*tr.stats.delta, tr.stats.delta)*1e3 #[ms]
    ii = np.where((st < t) & (t< et))[0]
    inds = np.intersect1d(range(BW, N-FW-1), ii)
#     print(t[inds])
    
    for i in inds:
#         print(i)
        FA = np.sum(d[range(i+1, i+FW+1, +1)])
        BA = np.sum(d[range(i-1, i-BW-1, -1)])
#         print(len(range(i+1, i+FW+1, +1)), len(range(i-1, i-BW-1, -1)))
        f[i] = (BW/FW) * (FA/BA) -1
    
    return f


def eval_biax_4m_eventloc(tpick, st, channel_loc, Vrng, PickNumStation=6, gridsize=1e-3, datacase="", vfixed_id=False):
    """
        evaluate location of event using grid search.
        
         [Input]
         stpickdat: Arrival time at each strain station (s)
        
         [Output]
         T: Occurence time (s)
         X: Hypocenter location X (m)
         Y: Hypocenter location Y (m)
         V: Velocity of wave speed (m/s)
         R: Error
        
        This script is transrated from [T,XY,V,R]=fmbiaxsthypolocGS(stpickdat) written by Futoshi Yamashita (2018/08/08).
        
        2021.4.7 Kurama Okubo
    """
    # Parameters fixed for 4m biax rock specimen
    L = 4.1 #[m] Length of rock
    W = 0.1 #[m] Width of rock
    
    # search stations used for relocation
    tpick_sorted = sorted(tpick.items(), key=lambda x: np.inf if np.isnan(x[1]) else x[1], reverse=False) # ignore the nan

    stlist = []
    xloc = np.zeros(0)
    yloc = np.zeros(0)
    zloc = np.zeros(0)
    pick = np.zeros(0)

    for i, (stnm, t1) in enumerate(tpick_sorted[:PickNumStation]):
        if np.isnan(t1):
            # skip this station
            continue;
                        
        stlist.append(stnm)
        xloc = np.append(xloc, channel_loc[stnm][0]/1e3) #[m]
        yloc = np.append(yloc, channel_loc[stnm][1]/1e3) #[m]
        zloc = np.append(zloc, channel_loc[stnm][2]/1e3) #[m]
        pick = np.append(pick, t1) #[ms]
    

    if len(stlist) < PickNumStation:
        warnings.warn("Number of stations {} for relocation is less than PickNumStation.".format(len(stlist)), UserWarning)
            
    # === Search range for XY === #
    X0 = xloc[0]
    # Xrng=np.linspace(X0-0.2, X0+0.2, int(np.round(0.4/gridsize)) + 1, endpoint=True) #[m] 
    Xrng=np.linspace(X0-0.1, X0+0.1, int(np.round(0.2/gridsize)) + 1, endpoint=True) #[m] 2022/11/7 update location 2021/12/22 narrow down area
    Xrng = np.delete(Xrng, np.where(Xrng<0.003));
    Xrng = np.delete(Xrng, np.where(Xrng>L-0.003));
    Nx=len(Xrng);

    Yrng=np.linspace(-W/2, W/2, int(np.round(W/gridsize)) + 1, endpoint=True) #[m] 
    Ny=len(Yrng);

    Nv=len(Vrng)

    Ermat=np.zeros((Nx,Ny));
    Vmat=np.zeros((Nx,Ny));
    Tmat=np.zeros((Nx,Ny));

    Erdat_all=np.zeros((Nx,Ny,Nv));

    # ===== Grid search =====
    for n1 in tqdm(range(Nx), desc="Grid searching:"):
        for n2 in range(Ny):
            Dmat=np.sqrt((xloc-Xrng[n1])**2+(yloc-Yrng[n2])**2+zloc**2);
            Erdat=np.ones(Nv);
            Tdat=np.ones(Nv);
            
            for n3 in range(Nv):
                Erdat[n3]=np.var(pick - 1e3*Dmat/Vrng[n3], ddof=1); #[ms]
                Tdat[n3]=np.mean(pick - 1e3*Dmat/Vrng[n3]); #[ms]

            Erdat_all[n1, n2, :] = Erdat

            if not vfixed_id:
                # grid searching v
                Ermat[n1,n2] = np.min(Erdat)
                tmpidx = np.argmin(Erdat)
                Vmat[n1,n2]=Vrng[tmpidx]
                Tmat[n1,n2]=Tdat[tmpidx]

            else:
                # fix v
                Ermat[n1,n2] = Erdat[vfixed_id]
                tmpidx = vfixed_id
                Vmat[n1,n2]=Vrng[tmpidx]
                Tmat[n1,n2]=Tdat[tmpidx]


    # ===== Optimal data =====
    
    [xid,yid]=np.unravel_index(np.argmin(Ermat), Ermat.shape)
    
    T = Tmat[xid,yid]
    X = Xrng[xid]
    Y = Yrng[yid]
    V=Vmat[xid,yid]
    R=Ermat[xid,yid]
    
    # debug
    if os.path.exists("../data/DEBUG_compute_CF"):
    	np.savez('../data/DEBUG_compute_CF/Erdat_{}.npz'.format(datacase), Erdat_all=Erdat_all, Xrng=Xrng, Yrng=Yrng, Vrng=Vrng)


    return (T, X, Y, V, R, Ermat, Vmat, Tmat, Xrng, Yrng)
