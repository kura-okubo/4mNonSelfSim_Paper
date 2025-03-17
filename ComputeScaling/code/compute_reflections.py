import numpy as np
from scipy.optimize import fsolve

def compute_p_reflection(tr_obs, param):
    # compute reflection from the side surfaces
    cp = param["cp"] #[m/s]
    rock_width = 0.1 #[m] 100mm width of fault
    rock_hight = 0.2 #[m] 200mm on onside of rock specimen
    sensorloc = np.array([0, 0, 70e-3])
    sourceloc = np.array([tr_obs.stats.xi1, tr_obs.stats.eta1, tr_obs.stats.zeta1])
    
    rvec = sourceloc - sensorloc
    # print(rvec)
    wvec_side = [0, -2*rock_width, 0]
    wvec_bottom = [0, 0, 2*(rock_hight-sensorloc[2])]
    wvec_top = [0, 0, -2*(rock_hight+sensorloc[2])]

    def f_lawcos(x):
        # compute law of cosines for the vector
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


def compute_ps_and_sp_side_reflection2(tr_obs, param):
    """
    compute ps and sp reflection from side surface by solving the simultaneous equation with Snell's law 
    """
    cp = param["cp"] # P wave velocity [m/s]
    cs = param["cs"] # S wave velocity[m/s]

    rock_width = 0.1 #[m] 100mm width of fault
    rock_hight = 0.2 #[m] 200mm on onside of rock specimen
    sensorloc = np.array([0, 0, 70e-3]) # location of sensor
    sourceloc = np.array([tr_obs.stats.xi1, tr_obs.stats.eta1, tr_obs.stats.zeta1]) # location of source

    rvec = np.array(sourceloc - sensorloc) # vector from sensor to source
    wvec_side = np.array([0, -2*rock_width, 0]) # vector from sensor to the virtual side of sensor for reflection

    # Function to compute pS
    vnorm = np.linalg.norm

    def func_res(x, type):
        """
        Define the simultaneous equations associated with the symmetry of reflection (isosceles triangle by mvec and wvec)
        and the Snell's law 
        """
        p, q = x # parameters to define the reflection point by mvec
        if p<0 or q<0:
            return [np.inf, np.inf] 

        mvec = p*rvec + q*wvec_side # vector from sensor to reflection point
        svec = mvec - rvec # vector from source to the reflection point
        cos_i = np.dot(svec, wvec_side)/(vnorm(svec) * vnorm(wvec_side)) # incident angle from source to vector 
        cos_j = np.dot(mvec, wvec_side)/(vnorm(mvec) * vnorm(wvec_side)) # reflected angle from source to vector 
        sin_i = np.sqrt(1-cos_i**2)
        sin_j = np.sqrt(1-cos_j**2)
        f = np.zeros(2)
        f[0] = np.dot(mvec, wvec_side) / vnorm(wvec_side)**2 - 0.5 # condition of isoscales triangle from sensor to the virtual side of sensor.
        if type=="pS":
            f[1] = sin_i/cp - sin_j/cs # Snell's law
        elif type=="sP":
            f[1] = sin_i/cs - sin_j/cp # Snell's law
        else:
            raise ValueError(f"type {type} is unknown.")
            
        return f

    tt_dict = dict()
    for type in ["pS", "sP"]:
        A = fsolve(func_res, [0.1, 0.1], args=(type,), full_output=True) # compute reflection point
        assert A[2] # check if the solution is found
        # print(A)
        p0, q0 = A[0]
        mvec_side = p0*rvec + q0*wvec_side # vector from sensor to the optimized reflection point
        svec_side = mvec_side - rvec # vector from source to the optimized reflection point
        ideg = np.rad2deg(np.arccos(np.dot(svec_side, wvec_side)/(vnorm(svec_side) * vnorm(wvec_side)))) # optimized incident angle
        jdeg = np.rad2deg(np.arccos(np.dot(mvec_side, wvec_side)/(vnorm(mvec_side) * vnorm(wvec_side)))) # optimized reflected angle

        l1 = np.sqrt( np.linalg.norm(rvec)**2 + np.linalg.norm(mvec_side)**2 - 2 * np.dot(rvec, mvec_side) ) # distance from source to reflection point
        l2 = np.sqrt( np.linalg.norm(mvec_side)**2 + np.linalg.norm(wvec_side)**2 - 2 * np.dot(mvec_side, wvec_side) ) # distance from reflection point to sensor
        
        if type=="pS":
            assert (ideg - jdeg) > 1e-8 # check if incident angle is larger than reflected angle
            tt_dict[type] = l1/cp + l2/cs # compute travel time as summation of incident and reflected parts

        elif type=="sP":
            assert (ideg - jdeg) < -1e-8 # check if incident angle is smaller than reflected angle
            tt_dict[type] = l1/cs + l2/cp # compute travel time as summation of incident and reflected parts

        else:
            raise ValueError(f"type {type} is unknown.")
    
        tt_dict[type+"_angles"] = (ideg, jdeg) # store incident and reflected angles

    return tt_dict["pS"], tt_dict["sP"], tt_dict["pS_angles"], tt_dict["sP_angles"]


# def compute_ps_and_sp_side_reflection(tr_obs, param):
#     """
#     compute ps and sp reflection from side surface by solving the simultaneous equation with Snell's law 
#     """
#     cp = param["cp"] #[m/s]
#     cs = param["cs"] #[m/s]

#     rock_width = 0.1 #[m] 100mm width of fault
#     rock_hight = 0.2 #[m] 200mm on onside of rock specimen
#     sensorloc = np.array([0, 0, 70e-3])
#     sourceloc = np.array([tr_obs.stats.xi1, tr_obs.stats.eta1, tr_obs.stats.zeta1])

#     rvec = np.array(sourceloc - sensorloc)
#     wvec_side = np.array([0, -2*rock_width, 0])

#     # Function to compute pS
#     # define m vector from sensor to reflection point
#     vnorm = np.linalg.norm

#     def func_ps(x):
#         p, q = x
#         if p<0 or q<0:
#             return [np.inf, np.inf] 
        
#         mvec = p*rvec + q*wvec_side
#         svec = mvec - rvec
#         cos_ps_i = np.dot(svec, wvec_side)/(vnorm(svec) * vnorm(wvec_side))
#         cos_ps_j = np.dot(mvec, wvec_side)/(vnorm(mvec) * vnorm(wvec_side))
#         sin_ps_i = np.sqrt(1-cos_ps_i**2)
#         sin_ps_j = np.sqrt(1-cos_ps_j**2)
#         f = np.zeros(2)
#         f[0] = np.dot(mvec, wvec_side) / vnorm(wvec_side)**2 - 0.5
#         f[1] = sin_ps_i/cp - sin_ps_j/cs
#         return f

#     A = fsolve(func_ps, [0.1, 0.1], full_output=True)
#     # print(A)
#     # compute i and j
#     p0_ps, q0_ps = A[0]
#     mvec_ps_side = p0_ps*rvec + q0_ps*wvec_side
#     svec_ps_side = mvec_ps_side - rvec

#     ideg_ps = np.rad2deg(np.arccos(np.dot(svec_ps_side, wvec_side)/(vnorm(svec_ps_side) * vnorm(wvec_side))))
#     jdeg_ps = np.rad2deg(np.arccos(np.dot(mvec_ps_side, wvec_side)/(vnorm(mvec_ps_side) * vnorm(wvec_side))))
#     assert (ideg_ps - jdeg_ps) > 1e-7

#     # compute distances
#     l1_ps = np.sqrt( np.linalg.norm(rvec)**2 + np.linalg.norm(mvec_ps_side)**2 - 2 * np.dot(rvec, mvec_ps_side) )
#     l2_ps = np.sqrt( np.linalg.norm(mvec_ps_side)**2 + np.linalg.norm(wvec_side)**2 - 2 * np.dot(mvec_ps_side, wvec_side) )
#     tps = l1_ps/cp + l2_ps/cs

#     # Function to compute sP
#     def func_sp(x):
#         p, q = x
#         if p<0 or q<0:
#             return [np.inf, np.inf] 
        
#         mvec = p*rvec + q*wvec_side
#         svec = mvec - rvec
#         cos_sp_i = np.dot(svec, wvec_side)/(vnorm(svec) * vnorm(wvec_side))
#         cos_sp_j = np.dot(mvec, wvec_side)/(vnorm(mvec) * vnorm(wvec_side))
#         sin_sp_i = np.sqrt(1-cos_sp_i**2)
#         sin_sp_j = np.sqrt(1-cos_sp_j**2)
#         f = np.zeros(2)
#         f[0] = np.dot(mvec, wvec_side) / vnorm(wvec_side)**2 - 0.5
#         f[1] = sin_sp_i/cs - sin_sp_j/cp
#         return f

#     B = fsolve(func_sp, [0.1, 0.1], full_output=True, xtol=1e-08)
#     # print(B)
#     # compute i and j
#     p0_sp, q0_sp = B[0]
#     mvec_sp_side = p0_sp*rvec + q0_sp*wvec_side
#     svec_sp_side = mvec_sp_side - rvec

#     ideg_sp = np.rad2deg(np.arccos(np.dot(svec_sp_side, wvec_side)/(vnorm(svec_sp_side) * vnorm(wvec_side))))
#     jdeg_sp = np.rad2deg(np.arccos(np.dot(mvec_sp_side, wvec_side)/(vnorm(mvec_sp_side) * vnorm(wvec_side))))

#     # print(ideg_sp, jdeg_sp)
#     assert (jdeg_sp - ideg_sp) > 1e-7

#     # compute distances
#     l1_sp = np.sqrt( np.linalg.norm(rvec)**2 + np.linalg.norm(mvec_sp_side)**2 - 2 * np.dot(rvec, mvec_sp_side) )
#     l2_sp = np.sqrt( np.linalg.norm(mvec_sp_side)**2 + np.linalg.norm(wvec_side)**2 - 2 * np.dot(mvec_sp_side, wvec_side) )
#     tsp = l1_sp/cs + l2_sp/cp

#     return tps, tsp, (ideg_ps, jdeg_ps, ideg_sp, jdeg_sp)

